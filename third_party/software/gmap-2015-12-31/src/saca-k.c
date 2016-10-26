#include "saca-k.h"
#include <stdio.h>
#include "bool.h"


// Author: Ge Nong,  Email: issng@mail.sysu.edu.cn
// Department of Computer Science, Sun Yat-sen University, 
// Guangzhou, China
// Date: December 24, 2012
//
// This is the demo source code for the algorithm SACA-K presented in this article:
// G. Nong, Practical Linear-Time O(1)-Workspace Suffix Sorting for Constant Alphabets, 
// ACM Transactions on Information Systems, Scheduled to Appear in July 2013.
// A draft for this article can be retrieved from http://code.google.com/p/ge-nong/.

#include <stdlib.h>

// set only the highest bit as 1, i.e. 1000...
const unsigned int EMPTY=((unsigned int)1)<<(sizeof(unsigned int)*8-1); 

// get s[i] at a certain level
#define chr(i) ((level==0)?((unsigned char *)s)[i]:((int *)s)[i])

void getBuckets(unsigned char *s, 
  unsigned int *bkt, unsigned int n,
  unsigned int K, bool end) { 
  unsigned int i, sum=0;
  
  // clear all buckets .
  for(i=0; i<K; i++) bkt[i]=0; 
  
  // compute the size of each bucket .
  for(i=0; i<n; i++) bkt[s[i]]++; 
  
  for(i=0; i<K; i++) { 
    sum+=bkt[i]; 
    bkt[i]=end ? sum-1 : sum-bkt[i]; 
  }
}

void putSuffix0(unsigned int *SA, 
  unsigned char *s, unsigned int *bkt, 
  unsigned int n, unsigned int K, int n1) {
  unsigned int i, j;

  // find the end of each bucket.
  getBuckets(s, bkt, n, K, true);

  // put the suffixes into their buckets.
  for(i=n1-1; i>0; i--) {
    j=SA[i]; SA[i]=0;
    SA[bkt[s[j]]--]=j;
  }
  SA[0]=n-1; // set the single sentinel suffix.
}

void induceSAl0(unsigned int *SA,
  unsigned char *s, unsigned int *bkt,
  unsigned int n, unsigned int K, bool suffix) {
  unsigned int i, j;

  // find the head of each bucket.
  getBuckets(s, bkt, n, K, false);

  bkt[0]++; // skip the virtual sentinel.
  for(i=0; i<n; i++)
    if(SA[i]>0) {
      j=SA[i]-1;
      if(s[j]>=s[j+1]) {
        SA[bkt[s[j]]]=j;
        bkt[s[j]]++;
        if(!suffix && i>0) SA[i]=0;
      }
    }
}

void induceSAs0(unsigned int *SA,
  unsigned char *s, unsigned int *bkt,
  unsigned int n, unsigned int K, bool suffix) {
  unsigned int i, j;

  // find the end of each bucket.
  getBuckets(s, bkt, n, K, true);

  for(i=n-1; i>0; i--)
    if(SA[i]>0) {
      j=SA[i]-1;
      if(s[j]<=s[j+1] && bkt[s[j]]<i) {
        SA[bkt[s[j]]]=j;
        bkt[s[j]]--;
        if(!suffix) SA[i]=0;
      }
    }
}

void putSubstr0(unsigned int *SA,
  unsigned char *s, unsigned int *bkt,
  unsigned int n, unsigned int K) {
  unsigned int i, cur_t, succ_t;

  // find the end of each bucket.
  getBuckets(s, bkt, n, K, true);

  // set each item in SA as empty.
  for(i=0; i<n; i++) SA[i]=0;

  succ_t=0; // s[n-2] must be L-type.
  for(i=n-2; i>0; i--) {
    cur_t=(s[i-1]<s[i] ||
           (s[i-1]==s[i] && succ_t==1)
          )?1:0;
    if(cur_t==0 && succ_t==1) SA[bkt[s[i]]--]=i;
    succ_t=cur_t;
  }

  // set the single sentinel LMS-substring.
  SA[0]=n-1;
}

void putSuffix1(int *SA, int *s, int n1) {
  int i, j, pos, cur, pre=-1;
  
  for(i=n1-1; i>0; i--) {
    j=SA[i]; SA[i]=EMPTY;
    cur=s[j];
    if(cur!=pre) {
      pre=cur; pos=cur;
    }
    SA[pos--]=j;
  }
}

void induceSAl1(int *SA, int *s, 
  int n, bool suffix) {
  int h, i, j, step=1;
  
  for(i=0; i<n; i+=step) {
    step=1; j=SA[i]-1;
    if(SA[i]<=0) continue;
    int c=s[j], c1=s[j+1];
    bool isL=c>=c1;
    if(!isL) continue;

    // s[j] is L-type.

    int d=SA[c];
    if(d>=0) {
      // SA[c] is borrowed by the left
      //   neighbor bucket.
      // shift-left the items in the
      //   left neighbor bucket.
      int foo, bar;
      foo=SA[c];
      for(h=c-1; SA[h]>=0||SA[h]==EMPTY; h--)
      { bar=SA[h]; SA[h]=foo; foo=bar; }
      SA[h]=foo;
      if(h<i) step=0;

      d=EMPTY;
    }

    if(d==EMPTY) { // SA[c] is empty.
      if(c<n-1 && SA[c+1]==EMPTY) {
        SA[c]=-1; // init the counter.
        SA[c+1]=j;
      }
      else        
        SA[c]=j; // a size-1 bucket.
    }
    else { // SA[c] is reused as a counter.
        int pos=c-d+1;
        if(pos>n-1 || SA[pos]!=EMPTY) {
          // we are running into the right
          //   neighbor bucket.
          // shift-left one step the items
          //   of bucket(SA, S, j).
          for(h=0; h<-d; h++)
            SA[c+h]=SA[c+h+1];
          pos--;
          if(c<i) step=0;
        }
        else
          SA[c]--;

        SA[pos]=j;
    }

    int c2;
    bool isL1=(j+1<n-1) && (c1>(c2=s[j+2]) || (c1==c2 && c1<i));  // is s[SA[i]] L-type?
    if((!suffix || !isL1) && i>0) {
      int i1=(step==0)?i-1:i;
      SA[i1]=EMPTY;
    }
  }

  // scan to shift-left the items in each bucket 
  //   with its head being reused as a counter.
  for(i=1; i<n; i++) {
    j=SA[i];
    if(j<0 && j!=EMPTY) { // is SA[i] a counter?
      for(h=0; h<-j; h++)
        SA[i+h]=SA[i+h+1];
      SA[i+h]=EMPTY;
    }
  }
}

void induceSAs1(int *SA, int *s, 
  int n, bool suffix) {
  int h, i, j, step=1;
  
  for(i=n-1; i>0; i-=step) {
    step=1; j=SA[i]-1;
    if(SA[i]<=0) continue;
    int c=s[j], c1=s[j+1];
    bool isS=(c<c1) || (c==c1 && c>i);
    if(!isS) continue;

    // s[j] is S-type

    int d=SA[c];
    if(d>=0) {
      // SA[c] is borrowed by the right
      //   neighbor bucket.
      // shift-right the items in the
      //   right neighbor bucket.
      int foo, bar;
      foo=SA[c];
      for(h=c+1; SA[h]>=0||SA[h]==EMPTY; h++)
      { bar=SA[h]; SA[h]=foo; foo=bar; }
      SA[h]=foo;
      if(h>i) step=0;

      d=EMPTY;
    }

    if(d==EMPTY) { // SA[c] is empty.
      if(SA[c-1]==EMPTY) {
        SA[c]=-1; // init the counter.
        SA[c-1]=j;
      }
      else
        SA[c]=j; // a size-1 bucket.
    }
    else { // SA[c] is reused as a counter.
        int pos=c+d-1;
        if(SA[pos]!=EMPTY) {
          // we are running into the left
          //   neighbor bucket.
          // shift-right one step the items 
          //   of bucket(SA, S, j).
          for(h=0; h<-d; h++)
            SA[c-h]=SA[c-h-1];
          pos++;
          if(c>i) step=0;
        }
        else
          SA[c]--;

        SA[pos]=j;
    }

    if(!suffix) {
      int i1=(step==0)?i+1:i;
      SA[i1]=EMPTY;
    }
  }

  // scan to shift-right the items in each bucket
  //   with its head being reused as a counter.
  if(!suffix)
    for(i=n-1; i>0; i--) {
      j=SA[i];
      if(j<0 && j!=EMPTY) { // is SA[i] a counter?
        for(h=0; h<-j; h++)
          SA[i-h]=SA[i-h-1];
        SA[i-h]=EMPTY;
      }
    }
}

void putSubstr1(int *SA, int *s, int n) {
  int h, i, j;

  for(i=0; i<n; i++) SA[i]=EMPTY;

  int c, c1, t, t1;
  c1=s[n-2];
  t1=0; 
  for(i=n-2; i>0; i--) {
    c=c1; t=t1; 
    c1=s[i-1];
    t1=c1<c || (c1==c && t);
    if(t && !t1) {
      if(SA[c]>=0) {
        // SA[c] is borrowed by the right
        //   neighbor bucket.
        // shift-right the items in the
        //   right neighbor bucket.
        int foo, bar;
        foo=SA[c];
        for(h=c+1; SA[h]>=0; h++)
        { bar=SA[h]; SA[h]=foo; foo=bar; }
        SA[h]=foo;

        SA[c]=EMPTY;
      }

      int d=SA[c];
      if(d==EMPTY) { // SA[c] is empty.
        if(SA[c-1]==EMPTY) {
          SA[c]=-1; // init the counter.
          SA[c-1]=i;
        }
        else
          SA[c]=i; // a size-1 bucket.
      }
      else { // SA[c] is reused as a counter
          int pos=c+d-1;
          if(SA[pos]!=EMPTY) {
            // we are running into the left
            //   neighbor bucket.
            // shift-right one step the items 
            //   of bucket(SA, S, i).
            for(h=0; h<-d; h++)
              SA[c-h]=SA[c-h-1];
            pos++;
          }
          else
            SA[c]--;

          SA[pos]=i;
      }
    }
  }

  // scan to shift-right the items in each bucket
  //   with its head being reused as a counter.
  for(i=n-1; i>0; i--) {
    j=SA[i];
    if(j<0 && j!=EMPTY) { // is SA[i] a counter?
      for(h=0; h<-j; h++)
        SA[i-h]=SA[i-h-1];
      SA[i-h]=EMPTY;
    }
  }

  // put the single sentinel LMS-substring.
  SA[0]=n-1;
}

unsigned int getLengthOfLMS(unsigned char *s, 
  unsigned int n, int level, unsigned int x) {
  if(x==n-1) return 1;  
  
  unsigned int dist, i=1;  
  while(1) {
    if(chr(x+i)<chr(x+i-1)) break;
    i++;
  }  
  while(1) {
    if(x+i>n-1 || chr(x+i)>chr(x+i-1)) break;
    if(x+i==n-1 || chr(x+i)<chr(x+i-1)) dist=i;
    i++;
  }
  
  return dist+1;
}

unsigned int nameSubstr(unsigned int *SA, 
  unsigned char *s, unsigned int *s1, unsigned int n, 
  unsigned int m, unsigned int n1, int level) {
  unsigned int i, j, cur_t, succ_t;

  // init the name array buffer
  for(i=n1; i<n; i++) SA[i]=EMPTY;

  // scan to compute the interim s1
  unsigned int name, name_ctr=0;
  unsigned int pre_pos, pre_len=0;
  for(i=0; i<n1; i++) {
    bool diff=false;
    unsigned int len, pos=SA[i];
    unsigned int d;

    len=getLengthOfLMS(s, n, level, pos);
    if(len!=pre_len) diff=true;
    else
      for(d=0; d<len; d++)
        if(pos+d==n-1 || pre_pos+d==n-1 ||
           chr(pos+d)!=chr(pre_pos+d)) {
          diff=true; break;
        }

    if(diff) {
      name=i; name_ctr++;
      SA[name]=1; // a new name.
      pre_pos=pos; pre_len=len;
    }
    else
      SA[name]++; // count this name.

    SA[n1+pos/2]=name;
  }

  // compact the interim s1 sparsely stored 
  //   in SA[n1, n-1] into SA[m-n1, m-1].
  for(i=n-1, j=m-1; i>=n1; i--)
    if(SA[i]!=EMPTY) SA[j--]=SA[i];

  // rename each S-type character of the
  //   interim s1 as the end of its bucket
  //   to produce the final s1.
  succ_t=1;
  for(i=n1-1; i>0; i--) {
    int ch=s1[i], ch1=s1[i-1];
    cur_t=(ch1< ch || (ch1==ch && succ_t==1))?1:0;
    if(cur_t==1) {
      s1[i-1]+=SA[s1[i-1]]-1;
    }
    succ_t=cur_t;
  }

  return name_ctr;
}

void getSAlms(unsigned int *SA, 
  unsigned char *s, 
  unsigned int *s1, unsigned int n, 
  unsigned int n1, int level ) {
  unsigned int i, j, cur_t, succ_t;

  j=n1-1; s1[j--]=n-1;
  succ_t=0; // s[n-2] must be L-type
  for(i=n-2; i>0; i--) {
    cur_t=(chr(i-1)<chr(i) ||
          (chr(i-1)==chr(i) && succ_t==1))?1:0;
    if(cur_t==0 && succ_t==1) s1[j--]=i;
    succ_t=cur_t;
  }

  for(i=0; i<n1; i++) SA[i]=s1[SA[i]];
  
  // init SA[n1..n-1]
  for(i=n1; i<n; i++) SA[i]=level?EMPTY:0; 
}


void SACA_K(unsigned char *s, unsigned int *SA,
  unsigned int n, unsigned int K,
  unsigned int m, int level) {
  unsigned int i;
  unsigned int *bkt=NULL;

  fprintf(stderr,"SACA_K called with n = %u, K = %u, level %d\n",n,K,level);

  // stage 1: reduce the problem by at least 1/2.

  if(level==0) {
    bkt=(unsigned int *)malloc(sizeof(int)*K);
    putSubstr0(SA, s, bkt, n, K);
    induceSAl0(SA, s, bkt, n, K, false);
    induceSAs0(SA, s, bkt, n, K, false);
  }
  else {
    putSubstr1((int *)SA, (int *)s,(int)n);
    induceSAl1((int *)SA, (int *)s, n ,false);
    induceSAs1((int *)SA, (int *)s, n, false);
  }

  // now, all the LMS-substrings are sorted and 
  //   stored sparsely in SA.

  // compact all the sorted substrings into
  //   the first n1 items of SA.
  // 2*n1 must be not larger than n.
  unsigned int n1=0;
  for(i=0; i<n; i++) 
    if((!level&&SA[i]>0) || (level&&((int *)SA)[i]>0))
      SA[n1++]=SA[i];

  unsigned int *SA1=SA, *s1=SA+m-n1;
  unsigned int name_ctr;
  name_ctr=nameSubstr(SA,s,s1,n,m,n1,level);

  // stage 2: solve the reduced problem.

  // recurse if names are not yet unique.
  if(name_ctr<n1)
    SACA_K((unsigned char *)s1, SA1, 
          n1, 0, m-n1, level+1);
  else // get the suffix array of s1 directly.
    for(i=0; i<n1; i++) SA1[s1[i]]=i;

  // stage 3: induce SA(S) from SA(S1).

  getSAlms(SA, s, s1, n, n1, level);
  if(level==0) {
    putSuffix0(SA, s, bkt, n, K, n1);
    induceSAl0(SA, s, bkt, n, K, true);
    induceSAs0(SA, s, bkt, n, K, true);
    free(bkt);
  }
  else {
    putSuffix1((int *)SA, (int *)s, n1);
    induceSAl1((int *)SA, (int *)s, n, true);
    induceSAs1((int *)SA, (int *)s, n, true);
  }
}
