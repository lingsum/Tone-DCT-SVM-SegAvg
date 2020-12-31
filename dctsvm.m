function dctsvm
% Tone recognition using DCT, SVM, and segment averaging
% Programmed using Octave 5.2.0 equipped with 
% LIBSVM (https://www.csie.ntu.edu.tw/~cjlin/libsvm/)

clc
clear all

% =========================================================
% Research results
% =========================================================
  % Musical equipment options: 'b':bellyra 'r':recorder 'p':pianica
  
  % Segment length options
  segment=[1 2 4 8 16 32 64 128 256 512];
  
  % SVM options (see LIBSVM documentary for further explanation)
  optx=['-g 32 -c 6400']
  %svmopt=['-t 0 -q ',optx]; % Linear kernel (optx: -c)
  %svmopt=['-t 1 -q ',optx]; % Polynomial kernel (optx: -d -g -c)
  svmopt=['-t 2 -q ',optx]; % Radial basis kernel (optx: -g -c)
  
  % Frame blocking 64
  for k=1:6
    rb64(k)=traintest('b',64,segment(k),svmopt);
    rr64(k)=traintest('r',64,segment(k),svmopt);
    rp64(k)=traintest('p',64,segment(k),svmopt);
  endfor

  % Frame blocking 128
  for k=1:7
    rb128(k)=traintest('b',128,segment(k),svmopt);
    rr128(k)=traintest('r',128,segment(k),svmopt);
    rp128(k)=traintest('p',128,segment(k),svmopt);
  endfor

  % Frame blocking 256
  for k=1:8
    rb256(k)=traintest('b',256,segment(k),svmopt);
    rr256(k)=traintest('r',256,segment(k),svmopt);
    rp256(k)=traintest('p',256,segment(k),svmopt);
  endfor

  % Frame blocking 512
  for k=1:9
    rb512(k)=traintest('b',512,segment(k),svmopt);
    rr512(k)=traintest('r',512,segment(k),svmopt);
    rp512(k)=traintest('p',512,segment(k),svmopt);
  endfor

  % Recognition results based on musical equipment
  z=9;
  rb64(z)=0;rb128(z)=0;rb256(z)=0;
  rr64(z)=0;rr128(z)=0;rr256(z)=0;
  rp64(z)=0;rp128(z)=0;rp256(z)=0;

  rbellyra=[rb64;rb128;rb256;rb512]
  rrecorder=[rr64;rr128;rr256;rr512]
  rpianica=[rp64;rp128;rp256;rp512]

endfunction

%==========================================================
% Internal function
%==========================================================
function outacc=traintest(museq,frame,seg,svmopt)
% Training and testing
% museq - musical equipment: 'b':bellyra 'r':recorder
%                            'p':pianica
% frame - frame length at frame blocking
% seg   - segment length at segment averaging

  % WAV path
  if museq=='b'
    wpath='D:\Data\Penelitian\alatmusik\bellyra\b';
  elseif museq=='r'
    wpath='D:\Data\Penelitian\alatmusik\recorder\r';
  else
    wpath='D:\Data\Penelitian\alatmusik\pianica\p';
  endif
    
  % Feature extraction of training data
  c1=xctraining(wpath,'c',frame,seg);
  c2=xctraining(wpath,'d',frame,seg);
  c3=xctraining(wpath,'e',frame,seg);
  c4=xctraining(wpath,'f',frame,seg);
  c5=xctraining(wpath,'g',frame,seg);
  c6=xctraining(wpath,'a',frame,seg);
  c7=xctraining(wpath,'b',frame,seg);
  c8=xctraining(wpath,'ch',frame,seg);

  % Training data - One vs All
  tr1=[c1;c2;c3;c4;c5;c6;c7;c8];
  tr2=[c2;c1;c3;c4;c5;c6;c7;c8];
  tr3=[c3;c1;c2;c4;c5;c6;c7;c8];
  tr4=[c4;c1;c2;c3;c5;c6;c7;c8];
  tr5=[c5;c1;c2;c3;c4;c6;c7;c8];
  tr6=[c6;c1;c2;c3;c4;c5;c7;c8];
  tr7=[c7;c1;c2;c3;c4;c5;c6;c8];
  tr8=[c8;c1;c2;c3;c4;c5;c6;c7];
  traindata={tr1,tr2,tr3,tr4,tr5,tr6,tr7,tr8};
  
  % Training label
  label1=ones(10,1);
  label2=2*ones(70,1);
  trainlabel=[label1;label2];
  
  % Feature extraction of testing data
  ts1=xctesting(wpath,'c',frame,seg);
  ts2=xctesting(wpath,'d',frame,seg);
  ts3=xctesting(wpath,'e',frame,seg);
  ts4=xctesting(wpath,'f',frame,seg);
  ts5=xctesting(wpath,'g',frame,seg);
  ts6=xctesting(wpath,'a',frame,seg);
  ts7=xctesting(wpath,'b',frame,seg);
  ts8=xctesting(wpath,'ch',frame,seg); 
 
  % Testing data 
  testdata=[ts1;ts2;ts3;ts4;ts5;ts6;ts7;ts8];
  
  % Testing label
  testlabel=1;
  
  % Input label
  x=ones(1,30);
  inlabel=[x 2*x 3*x 4*x 5*x 6*x 7*x 8*x];
  
  % SVM training and testing
  outacc=svmtraintest(trainlabel,traindata,...
         testlabel,testdata,inlabel,svmopt);
    
endfunction

% =======================================================
function xctr=xctraining(wpath,tone,frame,seg)
% Feature extraction for training data
  for k=1:10
    xctr(k,:)=xcdct([wpath tone num2str(k) '.wav'],frame,seg);
  endfor  
endfunction

% =======================================================
function xcts=xctesting(wpath,tone,frame,seg)
% Feature extraction for testing data
  for k=11:40
    xcts(k-10,:)=xcdct([wpath tone num2str(k) '.wav'],frame,seg);
  endfor
endfunction

% =======================================================
function xcout=xcdct(wavfile,frame,seg)
% Feature extraction
% wavfile : wav file input
% frame   : frame length at frame blocking
% seg     : segmen length at segment averaging
% xcout   : feature extraction output

  % Read wav file
  y0=audioread(wavfile);

  % Remove silence part
  y1=y0;
  b0=0.5*max(abs(y1));
  bx=find(y1>b0 | y1<-b0);
  y1(1:bx(1))=[];

  % Remove transition part
  y2=y1;
  b1=floor(0.3*5000);
  y2(1:b1)=[];

  % Frame blocking
  y3=y2(1:frame);

  % Normalization
  y3=y3/max(abs(y3));

  % Windowing
  h=hamming(frame);
  y4=y3.*h;

  % DCT
  y5=abs(dct(y4));
  y5(1)=0;

  % Logarithmic scaling
  y5=log10(y5+1);

  % Frame warping 
  y51=y5(1:frame/2); y52=y5((frame/2)+1:frame);
  y6=(y51+fliplr(y52))/2;

  % Segment averaging
  if seg~=1
    y6=reshape(y6,seg,[]);
    y7=mean(y6);
  else
    y7=y6;
  endif

  % Output format setting
  xcout=y7';
  
endfunction

% =======================================================
function outaccuracy=svmtraintest(trainlabel,traindata,...
  testlabel,testdata,inlabel,svmopt)
% SVM evaluation

  % SVM training
  for m=1:8
    svmModel{m}=svmtrain(trainlabel,traindata{m},svmopt);
  endfor

  % SVM testing 
  for m=1:240
    for n=1:8
      [plabel(n),~,decval(n)]=svmpredict(...
      testlabel,testdata(m,:),svmModel{n},'-q');
    endfor
    [~,outlabel(m)]=max(decval);
  endfor

  % Result
  numcorrect=length(find(abs(inlabel-outlabel)==0));
  outaccuracy=(numcorrect/length(outlabel))*100;
  
endfunction

% =======================================================


