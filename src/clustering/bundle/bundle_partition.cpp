#include "../../../include/clustering/bundle/bundle_partition.h"
#include "../../../include/clustering/bundle/bundle_text.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <numeric>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
using adaptive_bundle::require;
bool fits(int a,int b,int limit) { return a>0 && b>0 && a<=limit/b; }
}

AdaptivePlan makeAdaptivePlan(const std::vector<int>& dims,std::vector<AdaptiveEdge> edges,
    int order,int bundle_cap,double cutoff,bool allow_merges,AdaptiveLimits limits) {
    require(order>=1,"adaptive requires positive bundle order");
    require(bundle_cap>0,"the adaptive bundle dimension cap must be positive");
    require(std::isfinite(cutoff) && cutoff>=0,"adaptive merge cutoff must be finite and nonnegative");
    require(dims.size()<size_t(INT_MAX),"too many adaptive sites");
    const int n=int(dims.size());
    std::vector<int> parent(n),size=dims;
    std::iota(parent.begin(),parent.end(),0);
    std::vector<std::set<int>> neighbors(n);
    std::set<std::pair<int,int>> unique;
    for(int d:dims) require(d>0 && d<=bundle_cap,
        "an indivisible site exceeds bundle_max_dim");
    for(auto& e:edges) {
        if(e.a>e.b) std::swap(e.a,e.b);
        require(e.a>=0 && e.a<e.b && e.b<n,"invalid adaptive edge indices");
        require(std::isfinite(e.strength) && e.strength>=0,"invalid adaptive edge strength");
        require(unique.emplace(e.a,e.b).second,"duplicate adaptive edge");
        neighbors[e.a].insert(e.b);neighbors[e.b].insert(e.a);
    }
    std::sort(edges.begin(),edges.end(),[](const AdaptiveEdge& a,const AdaptiveEdge& b) {
        if(a.strength!=b.strength)return a.strength>b.strength;
        return std::make_pair(a.a,a.b)<std::make_pair(b.a,b.b);
    });
    auto root=[&](int i) {
        while(parent[i]!=i) { parent[i]=parent[parent[i]];i=parent[i]; }
        return i;
    };
    AdaptivePlan plan;
    for(const auto& e:edges) {
        if(!allow_merges)break;
        if(e.strength==0 || e.strength<cutoff)continue;
        ++plan.eligible;
        int a=root(e.a),b=root(e.b);
        if(a==b)continue;
        if(a>b)std::swap(a,b);
        if(!fits(size[a],size[b],bundle_cap)) { ++plan.blocked_bundle;continue; }
        // A cluster is `order` bundles, and its dimension follows from this cap
        // and the order rather than being capped in its own right, so nothing
        // here looks past the merged bundle (terminology 4절, 2026-09-18).
        const int merged=size[a]*size[b];
        std::set<int> next=neighbors[a];
        next.insert(neighbors[b].begin(),neighbors[b].end());next.erase(a);next.erase(b);
        parent[b]=a;size[a]=merged;size[b]=0;
        for(int c:next) {
            neighbors[c].erase(a);neighbors[c].erase(b);neighbors[c].insert(a);
        }
        neighbors[a]=std::move(next);neighbors[b].clear();++plan.merged;
    }
    std::vector<int> labels(n,-1);
    for(int i=0;i<n;++i) {
        const int a=root(i);
        if(labels[a]<0) {
            labels[a]=int(plan.bundles.size());plan.bundles.emplace_back();
            plan.dimensions.push_back(size[a]);
        }
        plan.bundles[labels[a]].push_back(i);
    }
    std::set<std::pair<int,int>> links;
    if(order>=2) for(const auto& e:edges) {
        int a=labels[root(e.a)],b=labels[root(e.b)];
        if(a>b)std::swap(a,b);
        if(a!=b)links.emplace(a,b);
    }
    plan.edges.assign(links.begin(),links.end());
    return plan;
}

AdaptivePlan makeAdaptiveAutoPlan(const std::vector<int>& dims,std::vector<AdaptiveEdge> edges,
    int order,int bundleCap,double cutoff,double duration,double budget,bool allow_merges,AdaptiveLimits limits) {
    require(order>=1,"adaptive requires positive bundle order");
    require(bundleCap>0,"the adaptive bundle dimension cap must be positive");
    require(std::isfinite(cutoff) && cutoff>=0,"adaptive merge cutoff must be finite and nonnegative");
    require(std::isfinite(duration) && duration>=0,"adaptive auto duration must be finite and nonnegative");
    require(std::isfinite(budget) && budget>=0,"adaptive auto connection budget must be finite and nonnegative");
    require(dims.size()<size_t(INT_MAX),"too many adaptive sites");
    for(int d:dims)require(d>0 && d<=bundleCap,"an indivisible site exceeds bundle_max_dim");
    const int n=int(dims.size());
    struct Link { long double strength=0; bool allowed=false; };
    std::vector<std::map<int,Link>> graph(n);
    std::vector<int> parent(n),sizes=dims;std::iota(parent.begin(),parent.end(),0);
    std::set<std::pair<int,int>> unique;
    for(auto& e:edges) {
        if(e.a>e.b)std::swap(e.a,e.b);
        require(e.a>=0 && e.a<e.b && e.b<n,"invalid adaptive edge indices");
        require(std::isfinite(e.strength) && e.strength>=0,"invalid adaptive edge strength");
        require(unique.emplace(e.a,e.b).second,"duplicate adaptive edge");
        // Keep zero-strength allowed links too: after merging they may make a
        // bundle-pair eligible, exactly as the actual-distance any-pair rule does.
        graph[e.a][e.b]=graph[e.b][e.a]=Link{e.strength,e.allowed};
    }
    // A link is omitted only because order 1 has no inter-bundle expansion or
    // because it is outside the allowed set. Cluster dimension no longer omits
    // anything: it follows from the bundle cap and the order.
    //
    // That is also why no merge is weighed against the budget any more. A merge
    // used to be able to push a neighbouring pair past the cluster cap and force
    // its link out, so each one was tried against a running total of forced
    // omissions. With one cap a merge can only take links out of the omitted set
    // - the merged pair's own, and a disallowed link absorbed by an allowed one -
    // so the total never rises. The check that remained fired only on rounding:
    // 130.40000000001963460... against ...4621..., one step of a long double, and
    // on that it refused a merge that was free (measured 2026-09-19, 300,000
    // random plans). The final audit below is what holds the budget.
    auto requiredOmission=[&](const Link& e) { return order==1 || !e.allowed; };
    auto root=[&](int a) { while(parent[a]!=a) {parent[a]=parent[parent[a]];a=parent[a];}return a; };
    std::sort(edges.begin(),edges.end(),[](const AdaptiveEdge& a,const AdaptiveEdge& b) {
        if(a.strength!=b.strength)return a.strength>b.strength;
        return std::make_pair(a.a,a.b)<std::make_pair(b.a,b.b);
    });
    AdaptivePlan plan;
    for(const auto& e:edges) {
        if(!allow_merges)break;
        if(!e.allowed || e.strength==0 || e.strength<cutoff)continue;
        ++plan.eligible;
        int a=root(e.a),b=root(e.b);if(a==b)continue;if(a>b)std::swap(a,b);
        if(!fits(sizes[a],sizes[b],bundleCap)) {++plan.blocked_bundle;continue;}
        const int merged=sizes[a]*sizes[b];
        std::map<int,Link> next;
        for(int g:{a,b})for(const auto& item:graph[g]) {
            const int c=item.first;const auto& link=item.second;
            if(c==a || c==b)continue;
            auto& updated=next[c];updated.strength+=link.strength;updated.allowed|=link.allowed;
        }
        parent[b]=a;sizes[a]=merged;sizes[b]=0;
        for(const auto& item:next) {
            auto& neighbor=graph[item.first];neighbor.erase(a);neighbor.erase(b);neighbor[a]=item.second;
        }
        graph[a]=std::move(next);graph[b].clear();++plan.merged;
    }
    std::vector<int> labels(n,-1);
    for(int i=0;i<n;++i) {
        const int g=root(i);
        if(labels[g]<0) {
            labels[g]=int(plan.bundles.size());plan.bundles.emplace_back();plan.dimensions.push_back(sizes[g]);
        }
        plan.bundles[labels[g]].push_back(i);
    }
    struct Candidate {int a,b;long double strength;};
    std::vector<Candidate> optional;
    for(int a=0;a<n;++a)if(parent[a]==a)for(const auto& item:graph[a]) {
        const int b=item.first;const auto& link=item.second;
        if(a>=b || link.strength==0)continue;
        ++plan.positive_pairs;plan.total_strength+=link.strength;
        if(requiredOmission(link)) {
            ++plan.omitted_pairs;plan.omitted_strength+=link.strength;
            if(!link.allowed) {++plan.radius_excluded_pairs;plan.radius_excluded_strength+=link.strength;}
        } else optional.push_back({labels[a],labels[b],link.strength});
    }
    if(plan.omitted_strength*duration>budget) {
        std::ostringstream message;
        message<<"adaptive auto connection budget exceeded by mandatory omissions: score="
               <<plan.omitted_strength*duration<<" budget="<<budget;
        if(order==1)
            message<<"; order 1 omits every remaining inter-bundle link: permit larger bundles/eligible merges, increase bundle order, or explicitly increase the budget";
        else
            message<<"; enlarge/remove rdip, lower rdipcut, adjust the bundle cap or bundle order, or explicitly increase the budget";
        throw std::runtime_error(message.str());
    }
    // For this fixed partition, dropping weakest pairs first minimizes the
    // number retained under the sum budget (not runtime or coherence error).
    std::sort(optional.begin(),optional.end(),[](const Candidate& a,const Candidate& b) {
        if(a.strength!=b.strength)return a.strength<b.strength;
        return std::make_pair(a.a,a.b)<std::make_pair(b.a,b.b);
    });
    bool retainRest=false;
    for(const auto& e:optional) {
        const long double proposed=plan.omitted_strength+e.strength;
        if(!retainRest && proposed*duration<=budget) {
            plan.omitted_strength=proposed;++plan.omitted_pairs;
        } else {
            retainRest=true;plan.edges.emplace_back(e.a,e.b);
            require(e.strength<=std::numeric_limits<double>::max(),"adaptive bundle coupling overflow");
            plan.selected_connections.push_back({e.a,e.b,double(e.strength)});
        }
    }
    require(std::isfinite(plan.total_strength*duration) && std::isfinite(plan.omitted_strength*duration),
            "adaptive auto coupling score overflow");
    require(plan.omitted_strength*duration<=budget,"adaptive auto final budget audit failed");
    std::sort(plan.edges.begin(),plan.edges.end());
    std::sort(plan.selected_connections.begin(),plan.selected_connections.end(),[](const AdaptiveEdge& a,const AdaptiveEdge& b) {
        return std::make_pair(a.a,a.b)<std::make_pair(b.a,b.b);
    });
    return plan;
}
