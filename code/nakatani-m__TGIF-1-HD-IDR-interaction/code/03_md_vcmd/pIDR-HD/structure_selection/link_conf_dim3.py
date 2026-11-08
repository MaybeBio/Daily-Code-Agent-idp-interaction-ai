import numpy as np
import subprocess as sub
import argparse
import random
import itertools as iter
import time
from tqdm import tqdm

def argarser():
    parser=argparse.ArgumentParser(description='')
    parser.add_argument('-s',help='',type=int)
    parser.add_argument('-n',help='',type=int)
    parser.add_argument('-inp',help='',type=str)
    #parser.add_argumet('-params',help='',type=str)
    args=parser.parse_args()
    return args



class reensemble:
    def __init__(self,stage,series):
        self.stage=stage
        self.series=series


    def read_inp(self,path):
        f=open(path,"r")
        #self.num_d1_state,self.num_d2_state=[int(i) for i in f.readline().split() ]
        #self.d1_rc_min,self.d1_rc_max=[ float(i) for i in f.readline().split()]
        self.d1_samp_vst_min,self.d1_samp_vst_max=[int(i) for i in f.readline().split()]
        #self.d2_rc_min,self.d2_rc_max=[ float(i) for i in f.readline().split()]
        self.d2_samp_vst_min,self.d2_samp_vst_max=[int(i) for i in f.readline().split()]
        #self.d1_bin,self.d2_bin=[float (i) for i in f.readline().split()]
        self.d3_samp_vst_min,self.d3_samp_vst_max=[int(i) for i in f.readline().split()]

        # return   self.d1_samp_state , self.d2_samp_state

    def read_params(self):
        f=open("../{0}/params_next.dat".format(self.stage),"r")
        f.readline()
        dim=int(f.readline())
        self.rc1_min_lists=[]
        self.rc1_max_lists=[]
        self.rc2_min_lists=[]
        self.rc2_max_lists=[]
        self.rc3_min_lists=[]
        self.rc3_max_lists=[]

        self.num_d1_state=int(f.readline())
        for i in range(self.num_d1_state):
            txt=f.readline().split()
            self.rc1_min_lists.append(float(txt[0]))
            self.rc1_max_lists.append(float(txt[1]))

        self.num_d2_state=int(f.readline())
        for i in range(self.num_d2_state):
            txt=f.readline().split()
            self.rc2_min_lists.append(float(txt[0]))
            self.rc2_max_lists.append(float(txt[1]))

        self.num_d3_state=int(f.readline())
        for i in range(self.num_d3_state):
            txt=f.readline().split()
            self.rc3_min_lists.append(float(txt[0]))
            self.rc3_max_lists.append(float(txt[1]))

        self.rc1_min_lists=self.rc1_min_lists[::-1]
        self.rc1_max_lists=self.rc1_max_lists[::-1]
        self.rc2_min_lists=self.rc2_min_lists[::-1]
        self.rc2_max_lists=self.rc2_max_lists[::-1]
        self.rc3_min_lists=self.rc3_min_lists[::-1]
        self.rc3_max_lists=self.rc3_max_lists[::-1]

    def static_vstate(self):
        def read_distance(stage,series):
            txt=sub.run(["tail -n1   ../"+ str(stage)+"/n"+str(series)+"/distance.dat "],encoding='utf-8',stdout=sub.PIPE,shell=True)
            nrc=txt.stdout.split()
            nrc1=float(nrc[1])
            nrc2=float(nrc[2])
            nrc3=float(nrc[3])
            return nrc1,nrc2,nrc3

        #rc1_min_lists=np.arange(self.d1_rc_min,self.d1_rc_max-self.d1_bin/2,self.d1_bin/2)[::-1]
        #rc2_min_lists=np.arange(self.d2_rc_min,self.d2_rc_max-self.d1_bin/2,self.d2_bin/2)[::-1]

        # print(rc1_min_lists)
        self.d_mesh=[ [ [ [] for d1 in range(self.num_d1_state) ] for d2 in range(self.num_d2_state)] for d3 in range(self.num_d3_state)]
        # d2_mesh=[ [] for i in range(len(self.state_num[0])) ]

        for  stage  in tqdm(range(1,self.stage+1)):
            for series in tqdm(range(1,self.series+1)):
                nrc1,nrc2,nrc3=read_distance(stage,series)
                d1_vstate,d2_vstate,d3_vstate=self.num_d1_state,self.num_d2_state,self.num_d3_state
               # print(d1_vstate)
               # print(d2_vstate)
               #print(d3_vstate)
                #vstate1,vstate2=cal_vstate(ncr1,ncr2)
                for rc1_min , rc1_max in zip (self.rc1_min_lists,self.rc1_max_lists):
                    if rc1_min<= nrc1 <= rc1_max :
                        break
                    d1_vstate-=1

                for rc2_min, rc2_max  in zip(self.rc2_min_lists, self.rc2_max_lists):
                    if rc2_min <= nrc2 <= rc2_max:
                        break
                    d2_vstate-=1

                for rc3_min, rc3_max  in zip(self.rc3_min_lists, self.rc3_max_lists):
                    if rc3_min <= nrc3 <= rc3_max:
                        break
                    d3_vstate-=1

                #self.d_mesh[d2_vstate-1][d1_vstate-1].append((stage,series))
                #print("OK")
                #print(d1_vstate)
                #print(d2_vstate)
                #print(d3_vstate)
                #print(self.d_mesh)
                self.d_mesh[d3_vstate-1][d2_vstate-1][d1_vstate-1].append(stage+(series*10**-4))
                #print("OK")


    def choice_fil(self):
        gro_fil_no=[]

        d1_samp_uni_list=np.random.randint(self.d1_samp_vst_min,self.d1_samp_vst_max+1,10)
        d2_samp_uni_list=np.random.randint(self.d2_samp_vst_min,self.d2_samp_vst_max+1,10)
        d3_samp_uni_list=np.random.randint(self.d3_samp_vst_min,self.d3_samp_vst_max+1,10)


        rand_samp_vst=list(iter.product(d1_samp_uni_list,d2_samp_uni_list,d3_samp_uni_list))

        """
        else:
            d_samp_vst=list(iter.product(d1_samp_vst_range,d2_samp_vst_range))
            random.shuffle(d_samp_vst)
        """
        # for  d2 in  range(self.d2_samp_state[0],self.d2_samp_state[1]+1):
        #  for d1 in range(self.d1_samp_state[0],self.d1_samp_state[1]+1):

        samp_vstate_list=[]
        for d1,d2,d3 in rand_samp_vst:
            gro_list=self.d_mesh[d3-1][d2-1][d1-1]
            if gro_list==[]:
                continue
            gro_fil_no.append(random.choice(gro_list))
            samp_vstate_list.append((d1,d2,d3))
            #print("OK")
        return  samp_vstate_list,gro_fil_no


    def make_n_folder(self):
        sub.call (["mkdir ../"+str(self.stage)+"/next_conf"],shell=True)
        for num in range(1,self.series+1):
            sub.call([ "mkdir ../"+str(self.stage) + "/next_conf/n"+str(num) ],shell=True)

    """
    def write_vstate(self,n,samp_list):
        fo=open("n{0}/state.dat".format(n),"w")

        fo.write("2")
    """
    """
    def link_gro(self,hit,gro_fil_no):

        for i in gro_fil_no:
            if hit == self.series:break
            hit+=1
            stage=i // (self.series*10)
            series=i % (self.series*10)
            print('cd {0}{1}/next_conf/n{2}/'.format(work_dir.stdout[:-18],self.stage,hit))
            print( 'ln -s   {0}{1}/n{2}/run.gro ./ '.format(work_dir.stdout[:-18],stage,series) )
        return hit
    """

def main():
    st=time.time()
    args=argarser()
    #print("check1")
    rem=reensemble(args.s,args.n)
    #print("check2")
    #d1_samp_state,d2_samp_state=rem.read_inp(args.inp)
    rem.read_inp(args.inp)
    #print("check3")
    rem.read_params()
    #print("check4")
    rem.make_n_folder()
    #print("check5")
    rem.static_vstate()
    #print("check")
    #print("#!/bin/bash")
    hit=0
    work_dir=sub.run(["pwd"],stdout=sub.PIPE ,encoding='utf-8', shell=True)
    while hit < args.n:
        d_samp_state,gro_fil_no=rem.choice_fil()
        for vst, i in zip(d_samp_state,gro_fil_no):
            #print(vst)
            if hit == args.n:
                break
            hit+=1
            stage= int(i)
            series= int(round((i - int(i))*10**4,0))
            print('cd {0}{1}/next_conf/n{2}/'.format(work_dir.stdout[:-18],args.s,hit))
            print( 'ln -s   {0}{1}/n{2}/run.gro ./ '.format(work_dir.stdout[:-18],stage,series) )
            print( 'echo 3 > state.dat')
            print( 'echo {0} >> state.dat '.format(vst[0]) )
            print( 'echo {0} >> state.dat '.format(vst[1]) )
            print( 'echo {0} >> state.dat '.format(vst[2]) )
        # hit=rem.link_gro(hit,gro_fil_no)
        # hit=len(gro_fil_no)

    end=time.time()
    print("echo "+str(end-st))
if __name__=="__main__":
    main()
