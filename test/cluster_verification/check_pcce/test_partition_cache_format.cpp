// Exercise the actual file codec under ASan/UBSan, without MPI initialization.
#include "../../../src/clustering/pcce_partition_cache.cpp"
#include <cassert>
#include <fstream>
#include <iterator>

int rank=0,nprocess=1,seed=19;
bool verbosity=false;

int main(int argc,char** argv){
    assert(argc==2);
    const char* path=argv[1];
    Partition_info p{2,3,6,0};
    Point centers[3]={{-0.0,1.0,2.0},{3.25,4.125,5.5},{6.1,7.2,8.3}};
    int assigned[6]={0,1,2,0,1,2},received[6]={};
    Point loaded[3]{};
    std::string key("identity\0with-binary-data",25),head=prefix(key),reason;
    writeCache(path,head,&p,centers,assigned);
    assert(readCache(path,head,&p,loaded,received,reason));
    assert(!std::memcmp(centers,loaded,sizeof(centers)));
    assert(!std::memcmp(assigned,received,sizeof(assigned)));
    std::ifstream in(path,std::ios::binary);
    std::string good((std::istreambuf_iterator<char>(in)),{});in.close();
    auto replace=[&](const std::string& data){std::ofstream out(path,std::ios::binary|std::ios::trunc);out.write(data.data(),data.size());};
    for(size_t i=0;i<good.size();++i){
        std::string bad=good;bad[i]^=1;replace(bad);
        assert(!readCache(path,head,&p,loaded,received,reason));
        replace(good.substr(0,i));
        assert(!readCache(path,head,&p,loaded,received,reason));
    }
    replace(good+"extra");assert(!readCache(path,head,&p,loaded,received,reason));
    replace(good);assert(!readCache(path,prefix(key+"changed"),&p,loaded,received,reason));
    assert(readCache(path,head,&p,loaded,received,reason));
    unlink(path);
    assert(!readCache(path,head,&p,loaded,received,reason) && reason=="missing");
    printf("PASS binary64 roundtrip and %zu corruption/truncation checks\n",2*good.size()+3);
}
