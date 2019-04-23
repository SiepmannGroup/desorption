/**@file Use MPI to farm out a large number of jobs to different processors.

   ssh / pdsh / pbsdsh based job distribution is sometimes
   limited by system settings. This program is supposed to
   get around the problem using MPI.

   @param Option 1:
   MPI-2 Standard supports dynamic processes with MPI_Comm_spawn()
   or MPI_Comm_spawn_multiple(), see:

   [1] MPI Standard (v2.2):
   http://www.mpi-forum.org/docs/mpi22-report/mpi22-report.htm

   [2] Open MPI documentation (v1.4):
   http://www.open-mpi.org/doc/v1.4/

   The current issue is that even after setting
   "ompi_non_mpi" to true, the calling process appears to
   hang, which is likely due to the Open MPI package being
   compiled without "--enable-mpi-thread-multiple".

   Add "-D__MPI_SPAWN__" to compilation, and run with:

   mpirun -np 1 ./jfmpi.exe

   @caution Implementation incomplete!

   @param Option 2:
   fork() is not supported by the MPI standard, so it is
   unlikely to be portable between multiple MPI
   implementations.

   Within Open MPI, fork() will almost certainly fail when
   using the OS-bypass networks (Myrinet, InfiniBand). The
   mechanisms that these networks use for high-speed/low
   latency are quite problematic with fork(). If using TCP
   or shared memory, fork() *might* work.

   Add "-D__FORK__" to compilation, and run with:

   mpirun -mca mpi_warn_on_fork 0 -bynode -np xxx ./jfmpi.exe
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <mpi.h>

using namespace std;

namespace mpi {
    static const int MASTER=0,BUFFER_SIZE=256,NUMBER_COMMAND=8192,NUMBER_ARGUMENT=8;
    int splitCommand(char *cmdArg,char *&cmd,char *arg[],const char *token);
    class MpiError;
    class MpiInstance;
}

using namespace mpi;

int mpi::splitCommand(char *cmdArg,char *&cmd,char *arg[NUMBER_ARGUMENT],const char *token=" ") {
    int iarg=0;
    cmd=strtok(cmdArg,token);
    do {
        arg[iarg]=strtok(NULL,token);
    } while (arg[iarg++]!=NULL);
    return iarg-1;
}

class mpi::MpiError {
    int id;
public:
    MpiError(int errorId=-1): id(errorId) {}
    ~MpiError() {}
    int getId() const {return id;}
    char *getMessage() const {return NULL;}
};

class mpi::MpiInstance {
    int universeSize;
    MPI_Comm serviceComm,workerComm;
    vector<MPI_Comm> childCommunicators;
    static const int TAG=1;
public:
    MpiInstance(int *argc,char ***argv,int nrankperjob=1) {
        int err,flag,*universeSizeP;
        if ((err=MPI_Init(argc,argv))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        if ((err=MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_UNIVERSE_SIZE,&universeSizeP,&flag))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        if (flag) {
            universeSize=*universeSizeP;
        } else {
            universeSize=getWorldSize();
        }
        // create service communicators
        int color=1,worldRank=getWorldRank();
        if ((worldRank%nrankperjob!=0) && (worldRank!=1)) {
            color=MPI_UNDEFINED;
        }
        if ((err=MPI_Comm_split(MPI_COMM_WORLD,color,worldRank,&serviceComm))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        // create worker communicators
        color=worldRank/nrankperjob;
        if (!worldRank) color=MPI_UNDEFINED;
        if ((err=MPI_Comm_split(MPI_COMM_WORLD,color,worldRank,&workerComm))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
#ifdef __DEBUG_MPI__
        cout << "[" << getWorldRank() << "/" << getServiceRank() << "/" << getWorkerRank() << "] worldSize: " << getWorldSize() << "; universeSize: " << universeSize << "; flag: " << flag << endl << endl;
#endif
    }
    ~MpiInstance() {
        int err=MPI_Finalize();
        if (MPI_SUCCESS!=err) exit(err);
    }
    MPI_Fint getWorkerCommFortran() {
        return MPI_Comm_c2f(workerComm);
    }
    int getRank(MPI_Comm comm) const {
        if (comm==MPI_COMM_NULL) return -1;
        int err,rank;
        if ((err=MPI_Comm_rank(comm,&rank))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        return rank;
    }
    int getServiceRank() const {
        return getRank(serviceComm);
    }
    int getWorkerRank() const {
        return getRank(workerComm);
    }
    int getWorldRank() const {
        return getRank(MPI_COMM_WORLD);
    }
    int getSize(MPI_Comm comm) const {
        int size;
        MPI_Comm_size(comm,&size);
        return size;
    }
    int getServiceSize() const {
        return getSize(serviceComm);
    }
    int getWorkerSize() const {
        return getSize(workerComm);
    }
    int getWorldSize() const {
        return getSize(MPI_COMM_WORLD);
    }
    int getUniverseSize() const {return universeSize;}
    int abort(int errorCode) const {
        int err=MPI_Abort(MPI_COMM_WORLD,errorCode);
        if (err!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        return err;
    }
    template<class T> int send(T *buf,int msglen,MPI_Datatype dtype,int dest) const {
        int err=MPI_Send(buf,msglen,dtype,dest,TAG,serviceComm);
#ifdef __DEBUG_MPI__
        cout << "[" << getWorldRank() << "/" << getServiceRank() << "/" << getWorkerRank() << "] send: " << *buf << " (len=" << msglen <<") to " << dest << ". Total size=" << sizeof buf << ", element size=" << sizeof *buf <<  endl << endl;
#endif
        if (err!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        return err;
    }
    int sendJob(char *buf,int dest) const {
        return send(buf,strlen(buf)+1,MPI_CHAR,dest);
    }
    int sendServiceRank() const {
        int serviceRank=getServiceRank();
        return send(&serviceRank,1,MPI_INTEGER,MASTER);
    }
    template<class T> int receive(T *buf,int msglen,MPI_Datatype dtype,int src) const {
        MPI_Status status;
        int err,count;
        if ((err=MPI_Recv(buf,msglen,dtype,src,TAG,serviceComm,&status))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        if ((err=MPI_Get_count(&status,dtype,&count))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
#ifdef __DEBUG_MPI__
        cout << "[" << getWorldRank() << "/" << getServiceRank() << "/" << getWorkerRank() << "] received: " << *buf << " (len=" << msglen << ", count=" << count << ") from " << src << ". Total size=" << sizeof buf << ", element size= " << sizeof *buf << endl << endl;
#endif
        return count;
    }
    int receiveJob(char *buf) const {
        return receive(buf,BUFFER_SIZE,MPI_CHAR,MASTER);
    }
    int receiveServiceRank(int *buf) const {
        return receive(buf,1,MPI_INTEGER,MPI_ANY_SOURCE);
    }
    int bcastJob(char *buf) const {
        int msglen=BUFFER_SIZE,err;
        if ((err=MPI_Bcast(buf,msglen,MPI_CHAR,MASTER,workerComm))!=MPI_SUCCESS) {
            throw MpiError(err);
        }
        return err;
    }
    bool isServiceRank() const {
        return !(serviceComm==MPI_COMM_NULL);
    }
#ifdef __MPI_SPAWN__
    void spawnMultiple(int nJob,char cmdArg[][BUFFER_SIZE],int numProcs[]) {
        int err,nArg,valuelen,flag;
        char *cmd[nJob],*arg[nJob][NUMBER_ARGUMENT],**arg1[nJob],value[200];
        MPI_Info info[nJob];
        MPI_Comm child;
#ifdef __DEBUG_MAIN__
        cout << "In spawnMultiple: nJob = " << nJob << endl << endl;
#endif
        for (int iJob=0;iJob<nJob;iJob++) {
            nArg=splitCommand(cmdArg[iJob],cmd[iJob],arg[iJob]);
            arg1[iJob]=MPI_ARGV_NULL;
            if (MPI_SUCCESS!=(err=MPI_Info_create(&info[iJob]))) {
                cerr << "create not success" << endl;
                throw MpiError(err);
            }
            if (MPI_SUCCESS!=(err=MPI_Info_set(info[iJob],"ompi_non_mpi","1"))) {
                cerr << "set not success" << endl;
                throw MpiError(err);
            }
        }
        childCommunicators.push_back(child);
        cout << "before spawn, err = " << err << endl;
        err=MPI_Comm_spawn_multiple(nJob,cmd,arg1,numProcs,info,getRank(),serviceComm,&child,MPI_ERRCODES_IGNORE);
        //if (MPI_SUCCESS!=(err=MPI_Comm_spawn_multiple(nJob,cmd,(char ***)arg,numProcs,info,getRank(),serviceComm,&child,MPI_ERRCODES_IGNORE))) {
        cout << "after spawn, err = " << err << endl;
        if (MPI_SUCCESS!=(err)) {
            throw MpiError(err);
        }
        return;
    }
#endif
};

//End of MPI related stuff

#ifdef __MPI_SPAWN__
int main(int argc,char *argv[]) {
    MpiInstance jobFarming(&argc,&argv);
    int numProcs[NUMBER_COMMAND];
    char command[BUFFER_SIZE];
    char commandArg[NUMBER_COMMAND][BUFFER_SIZE];

    if (jobFarming.getSize()!=1) {
        cerr << "more than 1 distributing processes started! exiting..." << endl;
        jobFarming.abort(-1);
    }

    int nFile=1,iJob=0,iNode=0;
    char **fileList;
    if (argc==1 || argc==2) {
        nFile=1;
        fileList=new char *[nFile];
        fileList[0]="command.lst";
    } else {
        nFile=argc-2;
        fileList=argv+2;
    }

#ifdef __DEBUG_MAIN__
    cout << nFile << " files to be read: " << fileList[0] << endl << endl;
#endif
    for (int iFile=0;iFile<nFile;iFile++) {
        ifstream commandList(fileList[iFile]);
        int nProcPerNode,nProcPerJob;
        commandList >> nProcPerNode >> nProcPerJob;
        commandList.getline(command,BUFFER_SIZE);
#ifdef __DEBUG_MAIN__
        cout << "# processors per node: "<<  nProcPerNode << ", # processors per job: " <<  nProcPerJob << endl << endl;
#endif
        do {
            commandList.getline(command,BUFFER_SIZE);
            numProcs[iJob]=nProcPerJob;
            strcpy(commandArg[iJob++],command);
            if (iJob * nProcPerJob > iNode * nProcPerNode) {
                if (iNode>jobFarming.getSize() && strlen(command)!=0) {
                    cerr << "required more nodes than available" << endl;
                    jobFarming.abort(-1);
                }
                iNode++;
            }
#ifdef __DEBUG_MAIN__
            cout << "command: " << command << "\niJob: " << iJob << ", iNode: " << iNode << "; size of command=" << sizeof command << ", size of a unit=" << sizeof *command << endl << endl;
#endif
        } while (strlen(command)!=0);
        commandList.close();
    }
    jobFarming.spawnMultiple(iJob-1,commandArg,numProcs);

    return 0;
}

#else
#include <sys/stat.h>
bool should_stop() {
    struct stat sts;
    if (stat("STOP",&sts) == 0) {
        return true;
    } else
        return false;
}

#include <string>
#include <set>
class JobManager {
    string *runningJobs;
    set<string> finishedJobs;
    fstream fdone;
public:
    JobManager(char *fname,int njobs): fdone(fname) {
        if (fdone.good()) {
            for (string line;getline(fdone,line);){
                finishedJobs.insert(line);
            }
        } else {
            fdone.open(fname,ios_base::out);
        }
        fdone.clear();
        runningJobs=new string[njobs];
        for (int i=0;i<njobs;i++) runningJobs[i]="";
    }
    ~JobManager() {
        fdone.close();
        delete[] runningJobs;
    }
    void addRunningJob(char *job,int rank) {
        runningJobs[rank-1]=job;
    }
    void markFinishedJob(int rank,ofstream& flog) {
        if (runningJobs[rank-1].empty()) return;
        flog << "Job " << runningJobs[rank-1] << " completed." << endl << endl;
        finishedJobs.insert(runningJobs[rank-1]);
        fdone << runningJobs[rank-1] << endl;
    }
    int isDone(char *job) const {
        return finishedJobs.count(string(job));
    }
};

#include <unistd.h>
#ifdef __FORK__
#include <sys/types.h>
#include <sys/wait.h>
//#include <string>
//#include <sstream>
//#include <iterator>
//#include <algorithm>
pid_t runCommand(char *cmd) {
    pid_t pId=fork();
    if (pId==0) { //Child
// !!execv() does not take container of strings...
// !!reverse to standard C library..
//         istringstream cmds(string(cmd));
//         vector<string> argvs;
//         copy(istream_iterator<string>(cmds),istream_iterator<string>(),back_inserter<vector<string> >(argvs));
        char *arg;
        char *argvs[strlen(cmd)];
        splitCommand(cmd,arg,&argvs[1]);
        argvs[0]=arg;
#ifdef __DEBUG_MAIN__
        cout << "process run " << argvs[0] << endl << endl;
#endif
        execv(argvs[0],argvs);
    } else if (pId<0) { //Failed
        cerr << "Failed to fork!" << endl;
    }
    return pId;
}
#else
#include <topmon_fortran_interface.h>
extern "C" {
    void sim_system_setup_mpi(int *,int *,int *);
    void topmon_main_monola(char *,int);
}
#endif

#include <time.h>
int main(int argc,char *argv[]) {
    MpiInstance jobFarming(&argc,&argv,atoi(argv[1]));
    char rootwd[BUFFER_SIZE],job[BUFFER_SIZE];

    if (!getcwd(rootwd,sizeof(rootwd))) {
        cerr <<  "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] Getting current working directory failed." << endl << endl;
    }

#ifdef __DEBUG_MAIN__
    cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] nrankperjob: " << atoi(argv[1]) << "; first job file #: " << atoi(argv[2]) << "; checkjob_frequency: " << atoi(argv[3]) << endl << endl;
#endif

    if (jobFarming.getWorldRank()==MASTER) {
        char *logFile="jobs/log.txt",*finishedJobList="jobs/done.txt",*file_prefix="jobs/run",*file_suffix=".txt",fname[BUFFER_SIZE];

        sprintf(fname,"%s/%s",rootwd,logFile);
        ofstream flog(fname);
        flog <<  "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] cwd: " << rootwd << endl << endl;

        sprintf(fname,"%s/%s",rootwd,finishedJobList);
        JobManager jm(fname,jobFarming.getServiceSize()-1);

        int file_start_no=atoi(argv[2]),slaveRank;
        sprintf(fname,"%s/%s%d%s",rootwd,file_prefix,file_start_no++,file_suffix);
        flog << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] Adding jobs from " << fname << endl << endl;
        ifstream jobs(fname);

        struct timespec checkjob_frequency;
        checkjob_frequency.tv_sec = atoi(argv[3]);
        checkjob_frequency.tv_nsec = 0;

    // If no more job files, then really no more jobs to launch
    bool no_more_jobs = false;
        while (1) {
      do {
        // read next job
        strcpy(job,"");
        if (should_stop()) {
          // test if user has instructed job to stop
          for (int i=1;i<jobFarming.getServiceSize();i++) {
        jobFarming.receiveServiceRank(&slaveRank);
        jm.markFinishedJob(slaveRank,flog);
        jobFarming.sendJob(job,slaveRank);
          }
          break;
        }
        if (jobs.eof()) {
          // open a new job list if the current one has reached EOF
          jobs.close();
          sprintf(fname,"%s/%s%d%s",rootwd,file_prefix,file_start_no++,file_suffix);
          jobs.open(fname);
          if (jobs.good()) flog << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] Adding jobs from " << fname << endl << endl;
          else no_more_jobs = true;
        }
        if (jobs.fail()) {
          // cannot open file
          file_start_no--;
          nanosleep(&checkjob_frequency,NULL);
        } else {
          jobs.getline(job,BUFFER_SIZE);
          if (jm.isDone(job)) strcpy(job,""); // test if the job is already finished
        }

        // if no more job, then escape while-loop, close file, and wait for service ranks to finish.
        if(no_more_jobs) {
          for (int i=1;i<jobFarming.getServiceSize();i++) {
        jobFarming.receiveServiceRank(&slaveRank);
        jm.markFinishedJob(slaveRank,flog);
        jobFarming.sendJob(job,slaveRank);
          }
          break;
        }
      } while (strlen(job)==0);

      fprintf(stdout,"Master rank %i testing job= %s with len= %i.\n",jobFarming.getWorldRank(),job,strlen(job));

      if (strlen(job)==0) break;
      jobFarming.receiveServiceRank(&slaveRank);
      jm.markFinishedJob(slaveRank,flog);
      jm.addRunningJob(job,slaveRank);
      flog << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] service rank " << slaveRank << " do: " << job << endl << endl;
      jobFarming.sendJob(job,slaveRank);
        }
        flog.close();
    fprintf(stdout,"Master rank %i just closed jobs file.\n",jobFarming.getWorldRank());
    } else {
        char pathname[BUFFER_SIZE];
#ifdef __FORK__
        int exitStatus;
        pid_t jobid,returnID;
#else
        int workerRank,workerCommFortran,workerSize;
#endif
        while (1) {
            if (jobFarming.isServiceRank()) {
                jobFarming.sendServiceRank();
        fprintf(stdout,"Service rank %i waiting for next job.\n",jobFarming.getWorldRank());
                jobFarming.receiveJob(job);
#ifdef __DEBUG_MAIN__
                cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] received: " << job << endl << endl;
#endif
            }
            jobFarming.bcastJob(job);
#ifdef __DEBUG_MAIN__
            cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] received broadcasted job: " << job << endl << endl;
#endif
            if (strlen(job)>0) {
#ifdef __FORK__
                jobid=runCommand(job);
#ifdef __DEBUG_FORK__
                cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] run: " << job << ", with jobid: " << jobid << endl << endl;
#endif
                returnID=waitpid(jobid,&exitStatus,0);
                if (returnID==0) {
                    cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] still running: " << jobid << endl << endl;
                } else if (returnID==jobid) {
                    if (WIFEXITED(exitStatus))
                        cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] child ended normally; jobid: " << jobid << endl << endl;
                    else if (WIFSIGNALED(exitStatus))
                        cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] child ended because of an uncaught signal; jobid: " << jobid << endl << endl;
                    else if (WIFSTOPPED(exitStatus))
                        cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] child has stopped; jobid: " << jobid << endl << endl;
                } else {
                    cout << "[" << jobFarming.getWorldRank() << "/" << jobFarming.getServiceRank() << "/" << jobFarming.getWorkerRank() << "] waitpid error; jobid: " << jobid << ", returnID: " << returnID << endl << endl;
                }
#else
                sprintf(pathname,"%s/%s",rootwd,job);
                chdir(pathname);
                workerRank=jobFarming.getWorkerRank();
                workerCommFortran=jobFarming.getWorkerCommFortran();
                workerSize=jobFarming.getWorkerSize();
                sim_system_setup_mpi(&workerSize,&workerRank,&workerCommFortran);
                strcat(pathname,"/topmon.inp");
                topmon_main_monola(pathname,strlen(pathname));
#endif
            } else {
         fprintf(stdout,"Service rank %i has non-positive job length.\n",jobFarming.getWorldRank());
         break;
            }
        }

   fprintf(stdout,"Service rank %i escaped while-loop.\n",jobFarming.getWorldRank());
    }

    fprintf(stdout,"Rank %i about to return 0.\n",jobFarming.getWorldRank());

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}
#endif
