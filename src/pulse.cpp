#include "../include/pulse.h"
#include "../include/memory.h"



/* Low-level functions ------------------------------------------------------*/

Pulse* Pulse_init(){
    Pulse* pulse = (Pulse*)allocArray1d(1, sizeof(Pulse));
    pulse->npulse = 0;
    pulse->pulseiter = false;
    pulse->pulsename[0] = '\0';
    pulse->sequence = NULL;
    return pulse;
}


void Pulse_setNpulse(Pulse* pulse, int npulse){
    pulse->npulse = npulse;
}

void Pulse_setPulsename(Pulse* pulse, char* pulsename){
    strcpy(pulse->pulsename, pulsename);
}

void Pulse_setPulseiter(Pulse* pulse, bool pulseiter){
    pulse->pulseiter = pulseiter;
}

void Pulse_setSequence_fromName(Pulse* pulse){

    if (pulse->sequence == NULL){
        printf("Error, pulse->sequence is not allocated! \n");
        exit(1);
    }

    int npulse = Pulse_getNpulse(pulse);
    bool made = false;

    if (strcmp(pulse->pulsename, "Ramsey") == 0 && npulse == 0){
        generateSequenceRamsey(&(pulse->sequence));
        made = true;
    }

    if (strcmp(pulse->pulsename, "HahnEcho") == 0 && npulse == 1){
        generateSequenceHE(&(pulse->sequence));
        made = true;
    }

    if (strcmp(pulse->pulsename, "CPMG") == 0 && npulse > 1){
        generateSequenceCPMG(&(pulse->sequence), pulse->npulse);
        made = true;
    }

    if (strcmp(pulse->pulsename, "Equal") == 0 && npulse > 1){
        generateSequenceEqual(&(pulse->sequence), pulse->npulse);
        made = true;
    }

    if (!made){
        printf("Error, pulse->pulsename is not matched! or you have to use sequence tag\n");
        exit(1);
    }

}

void Pulse_setSequence_fromInput(Pulse* pulse, double* seqinput){

    // e.g. Input file format :
    // if pulse->npulse = 4
    // input : seqinput : [0.2, 0.4, 0.6, 0.8] //gives pulse timing

    if (pulse->sequence == NULL){
        printf("Error, pulse->sequence is not allocated! \n");
        exit(1);
    }

    int npulse = Pulse_getNpulse(pulse);

    for (int i=0; i<(npulse+1); i++){

        double start = 0.0;
        double end = 1.0;

        if (i==0 && i!=npulse){
            start = 0.0;
            end = seqinput[i];    
        }
        else if (i!=0 && i==npulse){
            start = seqinput[i];
            end = 1.0;
        }
        else if (i==0 && i==npulse){
            start = 0.0;
            end = 1.0;
        }
        else{
            start = seqinput[i-1];
            end = seqinput[i];
        }

        (pulse->sequence)[i][0] = start; 
        (pulse->sequence)[i][1] = end;
        (pulse->sequence)[i][2] = end - start; 

        // Find 3-rd value
        bool findSameDiffIdx = false;
        for (int j=0; j<i; j++){
            if ((pulse->sequence)[j][2] == (pulse->sequence)[i][2]){
                (pulse->sequence)[i][3] = (double)j;
                if ((int)(pulse->sequence)[j][3] != (int)(pulse->sequence)[i][3]){
                    printf("Error, the SameDifferenceIndex is different! \n");
                }
                findSameDiffIdx = true;
                break;
            }
        }
        if (findSameDiffIdx){;}
        else{ (pulse->sequence)[i][3] = (double)i; }
    }
}

int Pulse_getNpulse(Pulse* pulse){
    return pulse->npulse;
}

char* Pulse_getPulsename(Pulse* pulse){
    return pulse->pulsename;
}

bool Pulse_getPulseiter(Pulse* pulse){
    return pulse->pulseiter;
}

double** Pulse_getSequence(Pulse* pulse){
    // data format
    //  Array format:
    //   sequence is started with 0.0 and end with 1.0
    //   sequence[npulse+1][4] (defulat : NULL)
    //   sequence[ipulse][0] = Fraction of Previous pulse
    //   sequence[ipulse][1] = Fraction of current pulse
    //   sequence[ipulse][2] = Difference of Fractions "[i][1] - [i][0]"
    //   sequence[ipulse][3] = The Index that have the same difference value
    //                         if there is no the same difference value 
    //                         then give the current index
    //                         ( This would reduce the calculational time cost )
    return pulse->sequence;
}

void Pulse_report(Pulse* pulse){

    int      npulse    = Pulse_getNpulse(pulse);
    char*    pulsename = Pulse_getPulsename(pulse);
    bool     pulseiter = Pulse_getPulseiter(pulse);
    double** sequence  = Pulse_getSequence(pulse);

    printLineSection();
    printTitle("Structure Pulse");

    printStructElementChar("pulsename", Pulse_getPulsename(pulse));
    printStructElementInt("npulse", Pulse_getNpulse(pulse));
    
    printStructElementBool("pulseiter", Pulse_getPulseiter(pulse));
    printf("%27s * ture  - pulse is applied to each qubit\n"," ");
    printf("%27s * false - pulse is applied to qubit-array\n"," ");
    
    if (sequence != NULL){
        printSubTitle("Sequence");
        printLine();
        printf("%6s               %8s %8s %8s %8s \n"," ","start","end","diff","index\n");
        for (int i=0; i<npulse+1; i++){
            printf("%6s sequence[%d] : %8.3lf %8.3lf %8.3lf %8.3lf \n"
            ," ",i,sequence[i][0],sequence[i][1],sequence[i][2],sequence[i][3]);
        }
        printLine();
    }
    printf("\n");
    printLineSection();
}

// allocation
void Pulse_allocSequence(Pulse* pulse){
    int npulse = Pulse_getNpulse(pulse);
    pulse->sequence = allocDouble2d(npulse+1, 4);
}

void Pulse_freeSequence(Pulse* pulse){
    int npulse = Pulse_getNpulse(pulse);
    freeDouble2d(pulse->sequence, npulse+1);
}

void Pulse_freeAll(Pulse* pulse){
    // Pulse_freeSequence(pulse);
    freeArray1d(pulse);
}


void generateSequenceRamsey(double*** sequence){

    (*sequence)[0][0] = 0.0; 
    (*sequence)[0][1] = 1.0;
    (*sequence)[0][2] = 1.0; 
    (*sequence)[0][3] = 0.0; 

}

void generateSequenceHE(double*** sequence){

    (*sequence)[0][0] = 0.0; 
    (*sequence)[0][1] = 0.5;
    (*sequence)[0][2] = 0.5; 
    (*sequence)[0][3] = 0.0; 

    (*sequence)[1][0] = 0.5; 
    (*sequence)[1][1] = 1.0;
    (*sequence)[1][2] = 0.5; 
    (*sequence)[1][3] = 0.0; 
   
}

void generateSequenceCPMG(double*** sequence, int npulse){

    // n : pulse
    // c : sequence i
    // sequence (2c-1)/2n * T

    double start = 0.0;
    double end = 1.0;

    for (int i=0; i<(npulse+1); i++){

        if (i==0 && i!=npulse){
            start = 0.0;
            end = double(2*(i+1)-1)/double(2*npulse); 
        }
        else if (i!=0 && i==npulse){
            start = double(2*i-1)/double(2*npulse);
            end = 1.0;
        }
        else if (i==0 && i==npulse){
            start = 0.0;
            end = 1.0;
        }
        else{
            start = double(2*i-1)/double(2*npulse);
            end = double(2*(i+1)-1)/double(2*npulse); 
        }

        (*sequence)[i][0] = start; 
        (*sequence)[i][1] = end;
        (*sequence)[i][2] = end - start; 

        // Find 3-rd value
        if (i==0 || i==npulse){
            (*sequence)[i][3] = 0.0; 
        }
        else{
            (*sequence)[i][3] = 1.0; 
        }
    }
}

void generateSequenceEqual(double*** sequence, int npulse){

    // n : pulse
    // c : sequence i
    // sequence c/n * T
    double start = 0.0;
    double end = 1.0;
    double diff = 1.0/double(npulse+1);

    for (int i=0; i<(npulse+1); i++){

        (*sequence)[i][0] = diff * i; 
        (*sequence)[i][1] = diff * (i+1);
        (*sequence)[i][2] = diff;
        (*sequence)[i][3] = 0.0; 

    }
}