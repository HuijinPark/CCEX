import time
import datetime
import psutil
import os
import sys
import argparse
import re
import copy

def _error(meg_error):
    print("\n    %s"%("_"*60));
    print(f"\n        ERROR !    {meg_error}\n\n")
    sys.exit(1)

class cal_form:

    def __init__(self,summary,mkdate):
        self.i_time  = time.time();
        self.summary = summary;
        self.mkdate  = mkdate;

    def _start(self,further_code_info=None):
        localtime = time.asctime(time.localtime(self.i_time))
        print("\n    %s"%("="*60))
        print("\n    <Summary> %s\n"%self.summary)
        print("\n    Latest update :  %s"%self.mkdate) # format Mon Apr 26 18:17:01 2021
        print("\n    Current time  : ",localtime); 
        print("\n    %s"%("_"*60));
        print("\n")
        if further_code_info!=None:
            code_info_copy = copy.deepcopy(further_code_info)
            try:
                del code_info_copy['summary']
            except:
                pass
            try:
                del code_info_copy['mkdate']
            except:
                pass
            for k,v in code_info_copy.items():
                print(f"\n    {k:13s} : {v}"); 
        print("\n    %s"%("_"*60));
        print("\n")

    def _update(self,further_code_info=None):
        localtime = time.asctime(time.localtime(self.i_time))
        print("\n    %s"%("="*60))
        print("\n    <Summary> %s\n"%self.summary)
        print("\n    Latest update :  %s"%self.mkdate) # format Mon Apr 26 18:17:01 2021
        print("\n    Current time  : ",localtime); 
        print("\n    %s"%("_"*60));
        print("\n")
        if further_code_info!=None:
            code_info_copy = copy.deepcopy(further_code_info)
            try:
                del code_info_copy['summary']
            except:
                pass
            try:
                del code_info_copy['mkdate']
            except:
                pass
            for k,v in code_info_copy.items():
                print(f"\n    {k:13s} : {v}"); 
        print("\n    %s"%("_"*60));
        print("\n")


    def _end(self,meg_error='non'):
        runtime = time.time() - self.i_time
        wall_clock_time = str(datetime.timedelta(seconds=runtime))#.split(".")
        print("\n    %s"%("_"*60));

        if meg_error == 'non':
            #
            #    Memory check
            #
            memory_usage_dict = dict(psutil.virtual_memory()._asdict())
            memory_usage_percent = memory_usage_dict['percent']
            print(f"\n    memory_usage_percent  : {memory_usage_percent: 9.3f} %")
            #
            #    current process RAM usage
            #
            pid = os.getpid()
            current_process = psutil.Process(pid)
            current_process_memory_usage_as_KB = current_process.memory_info()[0] / 2.**20
            print(f"\n    Current memory KB     : {current_process_memory_usage_as_KB: 9.3f} KB")
            #
            #
            print("\n    Wall clock time       : %s(s)"%(wall_clock_time))
            print("\n    %s"%("="*60))
            print("        JOB DONE    ")
            print("    %s"%("="*60))
        else:
            print(f"\n        ERROR !    {meg_error}\n\n")
            sys.exit(1)

    def is_argv_error(self,argv_dict,add_meg='none'):

        self.parser  = argparse.ArgumentParser(description = self.summary,\
                                               formatter_class=argparse.RawTextHelpFormatter,\
                                               epilog = add_meg);
        keyList = argv_dict.keys()
        varlist = []
        optlist = []
        
        for item in keyList:
            p = re.compile('[^=]+')
            res = p.findall(argv_dict[item])
            inf = re.findall(r'\w+',res[0])
            _type = self._type_is(inf[0],item)
            _var  = inf[-1]
            _help = res[1]
            _default = []
            _choices = []
            if 'key' in inf:
                optlist.append(_var)
                _default, _choices = self._option_is_in(res[1],item,_type)
#                When you need a set of choices
#                (self.parser).add_argument('--%s'%_var, type=_type, default = _default\
#                                            , choices = _choices, help = _help)
                (self.parser).add_argument('--%s'%_var, type=_type, default = _default\
                                            , help = _help)

#                print(_var, _type, _choices, _default)
            else:
                varlist.append(_var)
                (self.parser).add_argument(_var, type=_type\
                                            ,metavar = _var\
                                            ,help = _help)

        args = (self.parser).parse_args()
        print("\tInput argument : \n")
        for i in range(len(varlist)):
            print("\t [%d] %-7s : %s"%(i+1,varlist[i],vars(args)[varlist[i]]))
        print("\n")
        for opt in optlist:
            print("\t opt. --%-7s : %s"%(opt,vars(args)[opt]))

        return args;

    def _type_is(self,_i_type,item):
        if _i_type == 's':
            return str
        elif _i_type == 'f':
            return float
        elif _i_type == 'i':
            return int
        elif _i_type == 'b':
            return bool
        else:
            meg_error = f'input argument{item} error : type is wrong'
            self._end(meg_error)
    
    def _option_is_in(self,sentance,item,_type):
        _default = []
        _choices = []
        opt = re.search(r'\[[\S ]+\]',sentance)
        try:
            opt = re.findall(r'[^\/ \]\[]+',opt.group())
        except:
            opt = []
        for i in opt:
            if re.findall('([o|O][p|P][t|T])',i):
                pass;
            elif re.findall('{[a-zA-Z]+}',i):
                _default = re.search(r'[\w0-9.]+',i).group()
                _choices.append(_default)
            else:
                _choices.append(i)
        map(_type,_default)
        _choices = list(map(_type,_choices))

        return (_default, _choices) 
