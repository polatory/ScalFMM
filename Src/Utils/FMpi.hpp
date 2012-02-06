// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#ifndef FMPI_HPP
#define FMPI_HPP


#include <cstdio>

#include "FGlobal.hpp"
#include "FNoCopyable.hpp"
#include "FMath.hpp"



/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


#ifdef SCALFMM_USE_MPI
#include <mpi.h>
#else
struct MPI_Request{};
struct MPI_Status{ int MPI_SOURCE; };

typedef int MPI_Datatype;
typedef long long MPI_Offset;
typedef int MPI_Comm;
typedef int MPI_Group;
typedef int MPI_File;
typedef int MPI_Op;
typedef int MPI_Info;

MPI_Status* MPI_STATUSES_IGNORE = 0;
MPI_Status* MPI_STATUS_IGNORE = 0;


enum{
    MPI_SUCCESS,
    MPI_LONG_LONG,
    MPI_LONG,
    MPI_DOUBLE,
    MPI_FLOAT,
    MPI_INT,
    MPI_COMM_WORLD,
    MPI_BYTE,
    MPI_SUM,
    MPI_THREAD_MULTIPLE,
    MPI_ANY_SOURCE,
    MPI_MODE_RDONLY,
    MPI_INFO_NULL,
};

int MPI_Comm_rank( MPI_Comm comm, int *rank ){ return 0; }
int MPI_Comm_size( MPI_Comm comm, int *size ){ return 0; }

int MPI_Probe( int source, int tag, MPI_Comm comm, MPI_Status *status ){ return 0; }

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm){ return 0; }
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group){ return 0; }
int MPI_Comm_free(MPI_Comm *comm){ return 0; }
int MPI_Group_free(MPI_Group *group){ return 0; }
int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup){ return 0; }
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){ return 0; }
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided ){ return 0; }

int MPI_Finalize(){ return 0; }
int MPI_Abort(MPI_Comm comm, int errorcode){ return 0; }

int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm){ return 0; }

int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm){ return 0; }

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request){ return 0; }

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status){ return 0; }

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request){ return 0; }

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm){ return 0; }

int MPI_Alltoallv(void *sendbuf, int *sendcnts, int *sdispls,
                  MPI_Datatype sendtype, void *recvbuf, int *recvcnts,
                  int *rdispls, MPI_Datatype recvtype, MPI_Comm comm){ return 0; }

int MPI_Wait(MPI_Request *request, MPI_Status *status){ return 0; }

int MPI_Waitany(int count, MPI_Request array_of_requests[], int *index,
               MPI_Status *status){ return 0; }

int MPI_Waitall(int count, MPI_Request array_of_requests[],
               MPI_Status array_of_statuses[]){ return 0; }

int MPI_Waitsome(int incount, MPI_Request array_of_requests[],
                int *outcount, int array_of_indices[],
                MPI_Status array_of_statuses[]){ return 0; }

int MPI_File_open(MPI_Comm comm, char *filename, int amode,
                  MPI_Info info, MPI_File *fh){ return 0; }

int MPI_File_read(MPI_File mpi_fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status){ return 0; }

int MPI_File_get_position(MPI_File mpi_fh, MPI_Offset *offset){ return 0; }

int MPI_File_get_size(MPI_File mpi_fh, MPI_Offset *size){ return 0; }

int MPI_File_read_at(MPI_File mpi_fh, MPI_Offset offset, void *buf,
                    int count, MPI_Datatype datatype, MPI_Status *status){ return 0; }

int MPI_Get_count( MPI_Status *status,  MPI_Datatype datatype, int *count ){ return 0; }

int MPI_File_close(MPI_File *mpi_fh){ return 0; }

int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                int dest, int sendtag,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int source, int recvtag,
                MPI_Comm comm, MPI_Status *status){ return 0; }

int MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype,
                       int dest, int sendtag, int source, int recvtag,
                       MPI_Comm comm, MPI_Status *status){ return 0; }

int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root,
               MPI_Comm comm ){ return 0; }

#endif

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FMpi
* Please read the license
*
*/

class FMpi {
public:
////////////////////////////////////////////////////////
// MPI Flag
////////////////////////////////////////////////////////
    enum FMpiTag {
        // FMpiTreeBuilder
        TagExchangeIndexs,
        TagSplittedLeaf,
        TagExchangeNbLeafs,
        TagSandSettling,

        // FQuickSort
        TagQuickSort,

        // FMM
        TagFmmM2M,
        TagFmmL2L,
        TagFmmP2P,

        // Bitonic,
        TagBitonicMin,
        TagBitonicMax,
        TagBitonicMinMess,
        TagBitonicMaxMess,

        // Last defined tag
        TagLast,
    };

////////////////////////////////////////////////////////
// FComm to factorize MPI_Comm work
////////////////////////////////////////////////////////

    /** This class is used to put all the usual method
      * related mpi comm
      */
    class FComm : public FNoCopyable {
        int rank;   //< rank related to the comm
        int nbProc; //< nb proc in this group

        MPI_Comm communicator;  //< current mpi communicator
        MPI_Group group;        //< current mpi group


        // reset : get rank and nb proc from mpi
        void reset(){
            FMpi::Assert( MPI_Comm_rank(communicator,&rank),  __LINE__ );
            FMpi::Assert( MPI_Comm_size(communicator,&nbProc),  __LINE__ );
        }

    public:
        /** Constructor : dup the comm given in parameter */
        explicit FComm(MPI_Comm inCommunicator ) {
            FMpi::Assert( MPI_Comm_dup(inCommunicator, &communicator),  __LINE__ , "comm dup");
            FMpi::Assert( MPI_Comm_group(communicator, &group),  __LINE__ , "comm group");

            reset();
        }

        /** Free communicator and group */
        virtual ~FComm(){
            FMpi::Assert( MPI_Comm_free(&communicator),  __LINE__ );
            FMpi::Assert( MPI_Group_free(&group),  __LINE__ );
        }

        /** To get the mpi comm needed for communication */
        MPI_Comm getComm() const {
            return communicator;
        }

        /** The current rank */
        int processId() const {
            return rank;
        }

        /** The current number of procs in the group */
        int processCount() const {
            return nbProc;
        }

        ////////////////////////////////////////////////////////////
        // Split/Chunk functions
        ////////////////////////////////////////////////////////////

        /** Get a left index related to a size */
        template< class T >
        T getLeft(const T inSize)  const {
            const double step = (double(inSize) / double(processCount()));
            return T(FMath::Ceil(step * double(processId())));
        }

        /** Get a right index related to a size */
        template< class T >
        T getRight(const T inSize)  const {
            const double step = (double(inSize) / double(processCount()));
            const T res = T(FMath::Ceil(step * double(processId()+1)));
            if(res > inSize) return inSize;
            else return res;
        }

        /** Get a right index related to a size and another id */
        template< class T >
        T getOtherRight(const T inSize, const int other)  const {
            const double step = (double(inSize) / double(processCount()));
            const T res = T(FMath::Ceil(step * double(other+1)));
            if(res > inSize) return inSize;
            else return res;
        }

        /** Get a left index related to a size and another id */
        template< class T >
        T getOtherLeft(const T inSize, const int other) const {
            const double step = (double(inSize) / double(processCount()));
            return T(FMath::Ceil(step * double(other)));
        }

        /** Get a proc id from and index */
        template< class T >
        int getProc(const int position, const T inSize) const {
            const double step = (double(inSize) / processCount());
            return int(position/step);
        }

        ////////////////////////////////////////////////////////////
        // Mpi interface functions
        ////////////////////////////////////////////////////////////


        /** Reduce a value for proc == 0 */
        template< class T >
        T reduceSum(T data) const {
            T result(0);
            FMpi::Assert( MPI_Reduce( &data, &result, 1, FMpi::GetType(data), MPI_SUM, 0, communicator ), __LINE__);
            return result;
        }

        /** Reduce an average */
        template< class T >
        T reduceAverageAll(T data) const {
            T result[processCount()];
            FMpi::Assert( MPI_Allgather( &data, 1, FMpi::GetType(data), result, 1, FMpi::GetType(data), getComm()),  __LINE__ );

            T average = 0;
            for(int idxProc = 0 ; idxProc < processCount() ;++idxProc){
                average += result[idxProc] / processCount();
            }
            return average;
        }

        /** Change the group size */
        void groupReduce(const int from , const int to){
            int procsIdArray[to - from + 1];
            for(int idxProc = from ;idxProc <= to ; ++idxProc){
                procsIdArray[idxProc - from] = idxProc;
            }

            MPI_Group previousGroup = group;
            FMpi::Assert( MPI_Group_incl(previousGroup, to - from + 1 , procsIdArray, &group),  __LINE__ );

            MPI_Comm previousComm = communicator;
            FMpi::Assert( MPI_Comm_create(previousComm, group, &communicator),  __LINE__ );

            MPI_Comm_free(&previousComm);
            MPI_Group_free(&previousGroup);

            reset();
        }
    };

////////////////////////////////////////////////////////
// FMpi methods
////////////////////////////////////////////////////////

    /*
    We use init with thread because of an openmpi error:

    [fourmi062:15896] [[13237,0],1]-[[13237,1],1] mca_oob_tcp_msg_recv: readv failed: Connection reset by peer (104)
    [fourmi056:04597] [[13237,0],3]-[[13237,1],3] mca_oob_tcp_msg_recv: readv failed: Connection reset by peer (104)
    [fourmi053:08571] [[13237,0],5]-[[13237,1],5] mca_oob_tcp_msg_recv: readv failed: Connection reset by peer (104)

    Erreur pour le proc1
    [[13237,1],1][btl_openib_component.c:3227:handle_wc] from fourmi062 to: fourmi056 error polling LP CQ with status LOCAL LENGTH ERROR status number 1 for wr_id 7134664 opcode 0  vendor error 105 qp_idx 3
    Tous on la meme erreur le 2e 1 est remplacé par le rang.
    */
    FMpi(int inArgc, char **  inArgv ) : communicator(0) {
        int provided = 0;
        FMpi::Assert( MPI_Init_thread(&inArgc,&inArgv, MPI_THREAD_MULTIPLE, &provided), __LINE__);
        communicator = new FComm(MPI_COMM_WORLD);
    }

    /** Delete the communicator and call mpi finalize */
    ~FMpi(){
        delete communicator;
        MPI_Finalize();
    }

    /** Get the global communicator */
    const FComm& global() {
        return (*communicator);
    }

    ////////////////////////////////////////////////////////////
    // Mpi Types meta function
    ////////////////////////////////////////////////////////////

    static MPI_Datatype GetType(long long&){
        return MPI_LONG_LONG;
    }

    static MPI_Datatype GetType(long int&){
        return MPI_LONG;
    }

    static MPI_Datatype GetType(double&){
        return MPI_DOUBLE;
    }

    static MPI_Datatype GetType(float&){
        return MPI_FLOAT;
    }

    static MPI_Datatype GetType(int&){
        return MPI_INT;
    }

    ////////////////////////////////////////////////////////////
    // Mpi interface functions
    ////////////////////////////////////////////////////////////

    /** generic mpi assert function */
    static void Assert(const int test, const unsigned line, const char* const message = 0){
        if(test != MPI_SUCCESS){
            printf("[ERROR-QS] Test failled at line %d, result is %d", line, test);
            if(message) printf(", message: %s",message);
            printf("\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, int(line) );
        }
    }

    /** Compute a left index from data */
    template <class T>
    static T GetLeft(const T inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        return T(ceil(step * inIdProc));
    }

    /** Compute a right index from data */
    template <class T>
    static T GetRight(const T inSize, const int inIdProc, const int inNbProc) {
        const double step = (double(inSize) / inNbProc);
        const T res = T(ceil(step * (inIdProc+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    /** Compute a proc id from index & data */
    template <class T>
    static int GetProc(const T position, const T inSize, const int inNbProc) {
        const double step = double(inSize) / double(inNbProc);
        return int(double(position)/step);
    }

    /** assert if mpi error */
    static void MpiAssert(const int test, const unsigned line, const char* const message = 0){
        if(test != MPI_SUCCESS){
            printf("[ERROR] Test failled at line %d, result is %d", line, test);
            if(message) printf(", message: %s",message);
            printf("\n");
            fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, int(line) );
        }
    }

private:
    /** The original communicator */
    FComm* communicator;
};


#endif //FMPI_HPP


