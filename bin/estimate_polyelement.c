#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>

void error_allocation(char *format, ...) {
  va_list pa;
  int n;
  char *s, c;
  float f;
  fprintf(stderr,"erreur d'allocation : ");
  va_start(pa, format);
  while (*format != '\0') {
    if ( *format == '%' ) {
      switch (*++format) {
        case '%' : putc('%', stderr); break;
        case 'c' : /* affichage d'un caractère */
                   c = va_arg(pa, int);
                   putc(c, stderr);
                   break;
        case 'd' : /* affichage d'un entier */
                   n = va_arg(pa, int);
                   fprintf(stderr, "%d", n);
                   break;
        case 'f' : /* affichage d'un float */
                   f = va_arg(pa, double);    /* !!!!! */
                   fprintf(stderr,"%f", f);
                   break;
        case 's' : /* affichage d'une chaîne */
                   s = va_arg(pa, char *);
                   for ( ; *s != '\0'; s++ )
                     putc(*s,stderr);
                   break;
      } /* end switch */
    }
    else
      putc(*format,stderr);
    format++;
  }
  va_end(pa);

}

char* c_allocate(int (nb)){
    char* ptr = (char*)calloc(nb,sizeof(char));
        if(ptr==NULL)error_allocation("c_allocate");
        return ptr;
}

int* i_allocate(int (nb)){
    int* ptr = (int*)calloc(nb,sizeof(int));
        if(ptr==NULL)error_allocation("c_allocate");
        return ptr;
}


/*peremt de realouer un tableau de caractere taille total du nouveau tableau : nb*/
char* c_reallocate(char*ancienne, int (nb)){
    char* ptr = (char*)realloc(ancienne,nb*sizeof(char));
        if(ptr==NULL)error_allocation("c_reallocate");
        return ptr;
}


int CheckPattern(char* pattern, char*sequence , int lenpattern){
 int nbcommon=0;
 int cmt;
 for(cmt=0;cmt<lenpattern;cmt++){
  if(pattern[cmt]==sequence[cmt])nbcommon++;
 }
 return nbcommon;
} 

char* readline(FILE* read, char** lines, int* nballocline){
 int bal=1, cmt=0;
 char c;
 while(bal==1){
 c=fgetc(read);
 fprintf(stderr, "%c", c);
 if(c==EOF)bal=3;
 else if(c=='\n')bal=2;
 else{
  if(cmt>=(*nballocline)-1){
    (*lines)=c_reallocate((*lines),(*nballocline)+100);
    (*nballocline)+=100;
  }
  (*lines)[cmt]=c;
  cmt++;
 }
 } 
  (*lines)[cmt]='\0';
 return bal;
}

int researchpattern(char* seq,char* pattern, int lenpattern){
/*search first position of pattern*/ 
 int lenseq=strlen(seq),cmt, *posexpl;
 posexpl=i_allocate(lenseq) ;
 for(cmt=0;cmt<lenseq-lenpattern;cmt++){
  posexpl[cmt]=CheckPattern(pattern, &(seq[cmt]), lenpattern);
 }
 
 fprintf(stderr, "\n");
 free(posexpl);
 return 0;

}

int main(int argc, char *argv[])
{
    int Cmt, balEndFile, nballocline, lenpattern;
    char* fileinput;
    char* pattern, *lines;
    FILE *readseq;
    /*keep argument*/ 
    fprintf(stderr,"begin program : estimate size of \n");
    if(argc!=3){
    fprintf(stderr, "command line : executable fileinput pattern\n");
    exit(2);
    }
    fileinput=argv[1];
    pattern=argv[2];
    lenpattern=strlen(pattern);
    for(Cmt=0;Cmt<argc;Cmt++){
     fprintf(stderr,"%s\n",argv[Cmt]);
    }
    readseq= fopen(fileinput, "r"); // read mode
    if(readseq==NULL){
     fprintf(stderr, "can't read %s \n exit\n",fileinput);
     exit(2);
    }
    balEndFile=1;
   // allocate memory for lines :
   nballocline=10;
   lines=c_allocate(nballocline);
    while(balEndFile!=3){
     balEndFile=readline(readseq, &lines, &nballocline) ;
     researchpattern(lines,pattern, lenpattern);
     fprintf(stdout, "%s %d\n", lines,nballocline );
    }
    fclose(readseq) ;
    free(lines) ;   
    return 0;
}
