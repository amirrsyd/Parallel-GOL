#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

/***********************************************************
  Helper functions 
***********************************************************/

//For exiting on error condition
void die(int lineNo);

//For trackinng execution
long long wallClockTime();


/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue );

void freeSquareMatrix( char** );

void printSquareMatrix( char**, int size );


/***********************************************************
   World  related functions
***********************************************************/

#define ALIVE 'X' 
#define DEAD 'O'

char** readWorldFromFile( char* fname, int* size );

int countNeighbours(char** world, int row, int col);

void evolveWorld(char** curWorld, char** nextWorld, int size, int numTasks, int rank, int iteration);

char** extractRows(char** currWorld, int startRow, int numRows, int size);

/***********************************************************
   Simple circular linked list for match records
***********************************************************/

typedef struct MSTRUCT {
    int iteration, row, col, rotation;
    struct MSTRUCT *next;
} MATCH;


typedef struct {
    int nItem;
    MATCH* tail;
} MATCHLIST;

MATCHLIST* newList();

void deleteList( MATCHLIST*);

void insertEnd(MATCHLIST*, int, int, int, int);

void printList(MATCHLIST*);


/***********************************************************
   Search related functions
***********************************************************/

//Using the compass direction to indicate the rotation of pattern
#define N 0 //no rotation
#define E 1 //90 degree clockwise
#define S 2 //180 degree clockwise
#define W 3 //90 degree anti-clockwise

char** readPatternFromFile( char* fname, int* size );

void rotate90(char** current, char** rotated, int size);

void searchPatterns(char** world, int wSize, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list);

void searchSinglePattern(char** world, int wSize, int interation,
        char** pattern, int pSize, int rotation, MATCHLIST* list);



/***********************************************************
   Main function
***********************************************************/


int main( int argc, char** argv)
{
    int numtasks, rank, sendCount, recvCount, source;
    char **curW, **nextW, **temp, dummy[20];
    char **patterns[4];
    int dir, iterations, iter;
    int size, patternSize;
    long long before, after;
    MATCHLIST*list;

    if (argc < 4 ){
        fprintf(stderr, 
            "Usage: %s <world file> <Iterations> <pattern file>\n", argv[0]);
        exit(1);
    } 

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);



        curW = readWorldFromFile(argv[1], &size);
        nextW = allocateSquareMatrix(size+2, DEAD);

        char recvBuff[size];

        printf("World Size = %d\n", size);

        iterations = atoi(argv[2]);
        printf("Iterations = %d\n", iterations);



    patterns[N] = readPatternFromFile(argv[3], &patternSize);

    for (dir = E; dir <= W; dir++){
        patterns[dir] = allocateSquareMatrix(patternSize, DEAD);
        rotate90(patterns[dir-1], patterns[dir], patternSize);
    }
    printf("Rank = %d, ", rank);
    printf("Pattern size = %d\n", patternSize);
    printf("Rank: %d - After print pattern size\n", rank);
#ifdef DEBUG
    printSquareMatrix(patterns[N], patternSize);
    printSquareMatrix(patterns[E], patternSize);
    printSquareMatrix(patterns[S], patternSize);
    printSquareMatrix(patterns[W], patternSize);
#endif

 
    //Start timer
    before = wallClockTime();

    //Actual work start
    list = newList();

    /*
     *
     * ITERATION STARTS HERE!!!!
     *
     *
    */
    printf("Rank: %d - before iter...\n", rank);
    for (iter = 0; iter < iterations; iter++){
        if(rank==0){
#ifdef DEBUG
        printf("Rank: %d - World Iteration: %d\n", rank, iter);
        printSquareMatrix(curW, size+2);
#endif
            printf("Rank: %d - Searching pattern...\n", rank);
            searchPatterns( curW, size, iter, patterns, patternSize, list);
        }
        //Generate next generation
        printf("Rank: %d - What?!?! World is evolving!!!\n", rank);
        evolveWorld( curW, nextW, size, numtasks, rank, iter);
        temp = curW;
        curW = nextW;
        nextW = temp;
    }


    printList( list );

    //Stop timer
    after = wallClockTime();

    printf("Sequential SETL took %1.2f seconds\n", 
        ((float)(after - before))/1000000000);


    //Clean up
    deleteList( list );

    freeSquareMatrix( curW );
    freeSquareMatrix( nextW );

    freeSquareMatrix( patterns[0] );
    freeSquareMatrix( patterns[1] );
    freeSquareMatrix( patterns[2] );
    freeSquareMatrix( patterns[3] );

    MPI_Finalize();

    return 0;
}

/***********************************************************
  Helper functions 
***********************************************************/


void die(int lineNo)
{
    fprintf(stderr, "Error at line %d. Exiting\n", lineNo);
    exit(1);
}

long long wallClockTime( )
{
#ifdef __linux__
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}

/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue )
{

    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = (char*) malloc(sizeof(char) * size * size);
    if (contiguous == NULL) 
        die(__LINE__);


    memset(contiguous, defaultValue, size * size );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * size );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < size; i++){
        matrix[i] = &contiguous[i*size];
    }

    return matrix;
}

char** allocateRektMatrix( int size, int rows, char defaultValue )
{
    //printf("Rank: %d - Just entered rekt\n", rank);
    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = (char*) malloc(sizeof(char) * size * rows);
    if (contiguous == NULL) 
        die(__LINE__);


    memset(contiguous, defaultValue, size * rows );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * rows );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < rows; i++){
        matrix[i] = &contiguous[i*size];
    }

    //printf("Rank: %d - Leaving rekt\n", rank);
    return matrix;
}


void printSquareMatrix( char** matrix, int size )
{
    int i,j;
    
    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            printf("%c", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printRektMatrix(char** matrix, int size, int row){
    int i, j;
    for (i = 0; i < row; i++){
        for (j = 0; j < size; j++){
            printf("%c", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeSquareMatrix( char** matrix )
{
    if (matrix == NULL) return;

    free( matrix[0] );
}

/***********************************************************
   World  related functions
***********************************************************/

char** readWorldFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **world;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    //Using the "halo" approach
    // allocated additional top + bottom rows
    // and leftmost and rightmost rows to form a boundary
    // to simplify computation of cell along edges
    world = allocateSquareMatrix( size + 2, DEAD );

    for (i = 1; i <= size; i++){
        for (j = 1; j <= size; j++){
            fscanf(inf, "%c", &world[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }

    *sizePtr = size;    //return size
    return world;
    
}

int countNeighbours(char** world, int row, int col)
//Assume 1 <= row, col <= size, no check 
{
    int i, j, count;

    count = 0;
    for(i = row-1; i <= row+1; i++){
        for(j = col-1; j <= col+1; j++){
            count += (world[i][j] == ALIVE );
        }
    }

    //discount the center
    count -= (world[row][col] == ALIVE);

    return count;

}

void evolveWorld(char** curWorld, char** nextWorld, int size, int numTasks, int rank, int iteration)
{
    int i, j, q, liveNeighbours, numRows;
    char** smallNext;
    MPI_Status Stat;

    q = size/numTasks + (size % numTasks != 0);
    numRows = q + 2;
    int count = numRows*(size +2) ;

    char* recvBuff;
    recvBuff = (char*) malloc(sizeof(char*) * count);


    smallNext = allocateRektMatrix(size+2, numRows, DEAD);
    //recvBuff = allocateRektMatrix(size, numRows, DEAD);

    printf("Rank: %d - entering root conditional\n", rank);
    if(rank==0){
        for (i = 0; i < numTasks; i++){
            //toSend = extractRows(curWorld, i, numRows, size);
            printf("Rank: %d - root going to send in iteration %d\n", rank, iteration);
            //printRektMatrix(curWorld[q*1], size, numRows);
            MPI_Send(&curWorld[q*i][0], count, MPI_CHAR, i, (iteration*numTasks*2)+i, MPI_COMM_WORLD);
            printf("Rank: %d - root send in iteration %d liao\n", rank, iteration);
        }
    }

    printf("Rank: %d - Slave waiting in iteration %d\n", rank, iteration);
    MPI_Recv(&recvBuff[0], count, MPI_CHAR, 0, (iteration*numTasks*2)+rank, MPI_COMM_WORLD, &Stat);
    printf("Rank: %d - Slave received in iteration %d!\n", rank, iteration);

    char** workOn;

    workOn = allocateRektMatrix(size+2, numRows, DEAD);

    workOn[0] = recvBuff;
    for (i = 1; i < numRows; i++){
        workOn[i] = &recvBuff[i*(size+2)];
    }

    printf("Rank : %d - Subworld received in iteration %d:\n", rank, iteration);
    printRektMatrix(workOn, size+2, numRows);

    for (i = 1; i <= q; i++){
        for (j = 1; j <= size; j++){
            liveNeighbours = countNeighbours(workOn, i, j);
            smallNext[i][j] = DEAD;

            //Only take care of alive cases
            if (workOn[i][j] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    smallNext[i][j] = ALIVE;

            } else if (liveNeighbours == 3)
                    smallNext[i][j] = ALIVE;
        } 
    }

    int tag = (iteration * numTasks * 2) + numTasks +rank;
    int count2 = q*(size+2);

    printf("Rank: %d - Slave going to send back in iteration %d. Tag = %d\n", rank, iteration, tag);
    MPI_Send(&smallNext[1][0], count2, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
    printf("Rank: %d - Slave sent back data in iteration %d. Tag = %d\n", rank, iteration, tag);

    if(rank == 0){
        for (i = 0; i<numTasks; i++){
            printf("Rank: %d - Root waiting for return from process %d in iteration %d. Tag: %d\n", rank, i, iteration, (iteration * numTasks * 2) + numTasks +i);
            MPI_Recv(&nextWorld[i*q+1][0], count2, MPI_CHAR, i, (iteration * numTasks * 2) + numTasks +i, MPI_COMM_WORLD, &Stat);
            printf("Rank: %d - Root received return from process %d in iteration %d! Tag: %d\n", rank, i, iteration, (iteration * numTasks * 2) + numTasks +i);
            //for(j = 0; j < q; j++){
            //    memcpy(&nextWorld[i*q+j+1], &getBack[j], sizeof(nextWorld[i*q+j+1]));
            //}
        }
    }

}

char** extractRows(char** currWorld, int startRow, int numRows, int size){
    char** newMatrix;

    newMatrix = allocateRektMatrix(size, numRows, DEAD);

    for(int i = 0; i<numRows; i++){
        memcpy(&newMatrix[i], &currWorld[startRow++], sizeof(newMatrix[i]));
    }

    return newMatrix;
}

/***********************************************************
   Search related functions
***********************************************************/

char** readPatternFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **pattern;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    pattern = allocateSquareMatrix( size, DEAD );

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            fscanf(inf, "%c", &pattern[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }
    
    *sizePtr = size;    //return size
    return pattern;
}


void rotate90(char** current, char** rotated, int size)
{
    int i, j;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            rotated[j][size-i-1] = current[i][j];
        }
    }
}

void searchPatterns(char** world, int wSize, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list)
{
    int dir;

    for (dir = N; dir <= W; dir++){
        searchSinglePattern(world, wSize, iteration, 
                patterns[dir], pSize, dir, list);
    }

}

void searchSinglePattern(char** world, int wSize, int iteration,
        char** pattern, int pSize, int rotation, MATCHLIST* list)
{
    int wRow, wCol, pRow, pCol, match;


    for (wRow = 1; wRow <= (wSize-pSize+1); wRow++){
        for (wCol = 1; wCol <= (wSize-pSize+1); wCol++){
            match = 1;
#ifdef DEBUGMORE
            printf("S:(%d, %d)\n", wRow-1, wCol-1);
#endif
            for (pRow = 0; match && pRow < pSize; pRow++){
                for (pCol = 0; match && pCol < pSize; pCol++){
                    if(world[wRow+pRow][wCol+pCol] != pattern[pRow][pCol]){
#ifdef DEBUGMORE
                        printf("\tF:(%d, %d) %c != %c\n", pRow, pCol,
                            world[wRow+pRow][wCol+pCol], pattern[pRow][pCol]);
#endif
                        match = 0;    
                    }
                }
            }
            if (match){
                insertEnd(list, iteration, wRow-1, wCol-1, rotation);
#ifdef DEBUGMORE
printf("*** Row = %d, Col = %d\n", wRow-1, wCol-1);
#endif
            }
        }
    }
}

/***********************************************************
   Simple circular linked list for match records
***********************************************************/

MATCHLIST* newList()
{
    MATCHLIST* list;

    list = (MATCHLIST*) malloc(sizeof(MATCHLIST));
    if (list == NULL)
        die(__LINE__);

    list->nItem = 0;
    list->tail = NULL;

    return list;
}

void deleteList( MATCHLIST* list)
{
    MATCH *cur, *next;
    int i;
    //delete items first

    if (list->nItem != 0 ){
        cur = list->tail->next;
        next = cur->next;
        for( i = 0; i < list->nItem; i++, cur = next, next = next->next ) {
            free(cur); 
        }

    }
    free( list );
}

void insertEnd(MATCHLIST* list, 
        int iteration, int row, int col, int rotation)
{
    MATCH* newItem;

    newItem = (MATCH*) malloc(sizeof(MATCH));
    if (newItem == NULL)
        die(__LINE__);

    newItem->iteration = iteration;
    newItem->row = row;
    newItem->col = col;
    newItem->rotation = rotation;

    if (list->nItem == 0){
        newItem->next = newItem;
        list->tail = newItem;
    } else {
        newItem->next = list->tail->next;
        list->tail->next = newItem;
        list->tail = newItem;
    }

    (list->nItem)++;

}

void printList(MATCHLIST* list)
{
    int i;
    MATCH* cur;

    printf("List size = %d\n", list->nItem);    


    if (list->nItem == 0) return;

    cur = list->tail->next;
    for( i = 0; i < list->nItem; i++, cur=cur->next){
        printf("%d:%d:%d:%d\n", 
                cur->iteration, cur->row, cur->col, cur->rotation);
    }
}

