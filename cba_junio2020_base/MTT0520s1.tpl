//========================================================================

GLOBALS_SECTION

//========================================================================
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");

//========================================================================

TOP_OF_MAIN_SECTION

//========================================================================
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 

//========================================================================

DATA_SECTION

//========================================================================
 init_int ntime  
 init_int nedades
 init_int ntallas
 init_matrix mdatos(1,ntime,1,13)
 init_vector Tallas(1,ntallas)
 init_matrix Ctot(1,ntime,1,ntallas)
 init_matrix Ncru(1,ntime,1,ntallas)
 init_vector msex(1,ntallas)
 init_vector Wmed(1,ntallas)
 
 //!! ad_comm::change_datafile_name("ModtallasLL.ctl");
 init_number sigmaR
 init_vector dt(1,3)
 init_vector Par_bio(1,7)
 init_vector cv_Par_bio(1,7)
 init_int    minedad
 init_number bprior //Hiperestabilidad

  number log_Lr_prior
  number log_sr_prior
  number log_b_prior
  number log_beta_prior
  number log_Linf_prior
  number log_k_prior
   
  !! log_Linf_prior = log(Par_bio(1));
  !! log_k_prior = log(Par_bio(2));
  !! log_Lr_prior = log(Par_bio(3));
  !! log_sr_prior = log(Par_bio(4));
  !! log_beta_prior= log(Par_bio(5));
  !! log_b_prior = log(bprior);

 init_number L50prior
 init_number s1prior
 init_number s2prior
 init_int opt_sel2 // opcion domo en flota
 init_int opt_sel3 // opcion domo en crucero

 number log_L50prior
 number log_s1prior
 number log_s2prior

 !! log_L50prior = log(L50prior);
 !! log_s1prior = log(s1prior);
 !! log_s2prior = log(s2prior);

 init_int    nbloques1
 init_vector ybloques1(1,nbloques1)
 init_int    nbloques2
 init_vector ybloques2(1,nbloques2)
 init_int    nqbloques
 init_vector yqbloques(1,nqbloques)
 init_int    nqbloquesc
 init_vector yqbloquesc(1,nqbloquesc)
 init_int    nqbloquesmph
 init_vector yqbloquesmph(1,nqbloquesmph)

 init_number      prior_qf 
 init_number      prior_qc 
 init_number      prior_qmph 

 number      log_prior_qf 
 number      log_prior_qc 
 number      log_prior_qmph 

 !!log_prior_qf=log(prior_qf);
 !!log_prior_qc=log(prior_qc); 
 !!log_prior_qmph =log(prior_qmph);

 init_vector cv_q(1,3) 

 init_int    opt_qf
 init_int    opt_qc
 init_int    opt_qmph
 
 init_int    opt1_fase // selectividad flota
 init_int    opt2_fase // selectividad cruceros

 init_int    opt_Lr
 init_int    opt_sr
 init_int    opt_beta
 init_int    opt_vb1
 init_int    opt_vb2
 
 init_int    opt_Ro       //Opcion para estimar o dejar fijo log_Ro
 init_number log_priorRo  // Valor fijo para log_Ro, perfiles de verosimilitud 
 
 init_int    opt_F
 init_int    opt_devRt
 init_int    opt_devNo//Condicion inicial (Si no estima significa poblaciónen equilibrio)
 init_int    opt_bpow //hiperestabilidad

 init_int    npbr
 init_vector pbr(1,npbr)
 init_int ntime_sim

 init_number  festim
 init_number  oprec
 
//========================================================================

INITIALIZATION_SECTION

//========================================================================
  log_Rmed       8.8
  log_Lr         log_Lr_prior
  log_sr         log_sr_prior
  log_L50        log_L50prior 
  log_L50c       log_L50prior 
  log_sigma1     log_s1prior 
  log_sigma2     log_s2prior
  log_sigma1c    log_s1prior 
  log_sigma2c    log_s2prior
  log_b          log_b_prior 
  log_beta       log_beta_prior 
  log_k          log_k_prior
  log_Linf       log_Linf_prior
  log_qflo       log_prior_qf
  log_qcru       log_prior_qc
  log_qmph       log_prior_qmph

//========================================================================

PARAMETER_SECTION

//========================================================================
// selectividad paramétrica a la talla común
// init_bounded_vector log_L50f(1,nbloques1,-5,8,opt1_fase)  
 init_vector log_L50(1,nbloques1,opt1_fase)  
 init_vector log_sigma1(1,nbloques1,opt1_fase)
 init_vector log_sigma2(1,nbloques1,opt_sel2)

 init_vector log_L50c(1,nbloques2,opt2_fase)  
 init_vector log_sigma1c(1,nbloques2,opt2_fase)
 init_vector log_sigma2c(1,nbloques2,opt_sel3)

// parametros reclutamientos y mortalidades)
 init_number log_Rmed(opt_Ro)
 init_bounded_dev_vector log_desv_Rt(1,ntime,-10,10,opt_devRt)
 init_bounded_vector log_desv_No(1,nedades,-10,10,opt_devNo)
 init_bounded_vector log_F(1,ntime,-20,0.7,opt_F) // log  mortalidad por pesca por flota

// capturabilidades
 init_vector log_qflo(1,nqbloques,opt_qf)
 init_vector log_qcru(1,nqbloquesc,opt_qc)
 init_vector log_qmph(1,nqbloquesmph,opt_qmph)

 init_number log_b(opt_bpow)

// Crecimiento
 init_number log_Lr(opt_Lr)
 init_number log_sr(opt_sr)
 init_number log_beta(opt_beta)
 init_number log_Linf(opt_vb1)
 init_number log_k(opt_vb2)

//---------------------------------------------------------------------------------
//Defino las variables de estado 
 vector BMflo(1,ntime)
 vector BMcru(1,ntime)
 vector Brec(1,ntime)
 vector BMmph(1,ntime)
 sdreport_vector pred_CPUE(1,ntime);
 sdreport_vector pred_Bcru(1,ntime);
 sdreport_vector pred_Desemb(1,ntime);
 vector likeval(1,9);
 vector Neq(1,ntallas);

 sdreport_vector Rpred(1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_year(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector delta(1,ntallas)
 vector Lesp(1,ntallas)
 vector sigmaL(1,ntallas)
 vector pre(1,ntallas)

 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 vector BDo(1,ntime);
 vector No(1,ntallas)
 vector prior(1,7)
 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUE(1,ntime);
 vector Bcru(1,ntime);
 vector dt_C(1,ntime);
 vector Frms_bloque(1,ntime);
 vector mph(1,ntime);
 vector Lmed_obs(1,ntime)
 sdreport_vector Lmed_pred(1,ntime)
 vector Lmed_obsc(1,ntime)
 sdreport_vector Lmed_predc(1,ntime)
 vector edades(1,nedades)
 sdreport_vector Reclutas(1,ntime)
 vector nm(1,ntime)
 vector nmc(1,ntime)
 vector penalty(1,7)


 matrix cv_index(1,4,1,ntime)

 matrix S1(1,nbloques1,1,ntallas)
 matrix S2(1,nbloques2,1,ntallas)
 matrix Sel(1,ntime,1,ntallas)
 matrix Selc(1,ntime,1,ntallas)
 matrix F(1,ntime,1,ntallas)
 matrix Z(1,ntime,1,ntallas)
 matrix S(1,ntime,1,ntallas)
 matrix N(1,ntime,1,ntallas)

 matrix NM(1,ntime,1,ntallas)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)
 matrix NVcru(1,ntime,1,ntallas)
 matrix NVmph(1,ntime,1,ntallas)
 matrix TEMP(1,ntime,1,ntallas)
 matrix pred_Ctot(1,ntime,1,ntallas)
 matrix pobs(1,ntime,1,ntallas)
 matrix ppred(1,ntime,1,ntallas)
 matrix pobsc(1,ntime,1,ntallas)
 matrix ppredc(1,ntime,1,ntallas)
 matrix T(1,ntallas,1,ntallas)
 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,ntallas)
 number suma1
 number suma2
 number suma3
 number suma4
 number suma5
  
 number So
 number alfa
 number beta

 number Linf
 number k
 number Linfh
 number M
 number Lr
 number sr
 number Lm
 number Rm
 number h

 number BDp
 number Npplus
 number Bp_anch 

 number nm1;
 number cuenta1;
 number alfa_sr;
 number beta_sr;
 number pF

 vector Npact(1,ntallas)
 vector Np(1,ntallas)
 vector Zpbr(1,ntallas)
 sdreport_vector Fpbr(1,ntallas)
 vector Sp(1,ntallas)

 sdreport_matrix Bp(1,npbr,1,ntime_sim)
 
 //cba actual
 vector NV(1,ntallas)
 vector CTP(1,ntallas)
 sdreport_vector Yact(1,npbr)
 //cba proyectada
 vector NVp(1,ntallas) 
 vector CTPp(1,ntallas)
 sdreport_matrix Yp(1,npbr,1,ntime_sim)
 
 matrix Rpp(1,npbr,1,ntime_sim)
 
 objective_function_value f
  
 sdreport_vector BD(1,ntime) // 
 sdreport_vector BT(1,ntime) // 
 sdreport_vector RPRlp(1,ntime) // 
 sdreport_vector RPR(1,ntime)
 sdreport_vector Frpr(1,ntime)
 sdreport_vector pred_mph(1,ntime);
 sdreport_number SSBo
 


  

  
//========================================================================

PRELIMINARY_CALCS_SECTION

//========================================================================
 yrs=column(mdatos,1);
 Desemb=column(mdatos,2);
 CPUE=column(mdatos,4);
 Bcru=column(mdatos,6);
 nm=column(mdatos,8);
 nmc=column(mdatos,9);
 mph=column(mdatos,10);
 dt_C=column(mdatos,12);
 Frms_bloque=column(mdatos,13);
 //cout << "dt_C:" << dt_C << endl;exit(1);

 edades.fill_seqadd(minedad,1);

 cv_index(1)=column(mdatos,3);
 cv_index(2)=column(mdatos,5);
 cv_index(3)=column(mdatos,7);
 cv_index(4)=column(mdatos,11);

 Linf=Par_bio(1);
 k=Par_bio(2);
 M=Par_bio(6);
 h=Par_bio(7);

 Unos_tallas=1;// lo uso en operaciones matriciales con tallas
 Unos_year=1;// lo uso en operaciones matriciales con el año

//========================================================================

RUNTIME_SECTION

//========================================================================
  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-3,1e-5
  maximum_function_evaluations 100,100,200,3000,3500

//========================================================================

PROCEDURE_SECTION

//========================================================================
// se listan las funciones que contienen los calculos
 Eval_Trans_talla_talla();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_logverosim();
 Eval_funcion_objetivo();

 if(last_phase()){Eval_CTP();}


//========================================================================

FUNCTION Eval_Trans_talla_talla

//========================================================================
  Linf=exp(log_Linf);
  k=exp(log_k);
  beta=exp(log_beta);

//  if(active(log_k)){k=mfexp(log_k);}
  if(active(log_beta)){beta=mfexp(log_beta);}

 int i, j;
 
// matriz de transicion modelo normal

  delta=(Linf-Tallas)*(1-mfexp(-k));// incremento en tallas
  Lesp=Tallas+delta; // talla esperada luego del crecimiento
  sigmaL=delta*beta;  

  for (i=1;i<=ntallas;i++){
    for (j=1;j<=ntallas;j++){
      if(i==j){
         T(i,j)=1.0;}}
   }


  for (i=1;i<=ntallas;i++){

    for (j=1;j<=ntallas;j++){
     if(sigmaL(i)>0){
     T(i,j)=mfexp(-0.5*square((Lesp(i)-Tallas(j))/sigmaL(i)));}}
   }


  for (j=1;j<=ntallas;j++){
  T(j)/=sum(T(j));
  } 

//========================================================================

FUNCTION Eval_selectividad

//========================================================================
 int i,j;

 // FLOTA...................

 for (j=1;j<=nbloques1;j++){

 S1(j)=exp(-0.5*square(Tallas-exp(log_L50(j)))/square(exp(log_sigma1(j))));


    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50(j))){
      S1(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50(j)))/square(exp(log_sigma2(j))));
      }

 }}

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sel(i)=S1(j);}
       }
   }

 // CRUCERO...................

 for (j=1;j<=nbloques2;j++){

 S2(j)=exp(-0.5*square(Tallas-exp(log_L50c(j)))/square(exp(log_sigma1c(j))));

    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50c(j))){
      S2(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50c(j)))/square(exp(log_sigma2c(j))));
      }

 }}


   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques2;j++){
              if (yrs(i)>=ybloques2(j)){
                Selc(i)=S2(j);}
       }
   }

//========================================================================

FUNCTION Eval_mortalidades

//========================================================================

 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_tallas));
 Z=F+M;
 S=mfexp(-1.0*Z);

//========================================================================

FUNCTION Eval_abundancia

//========================================================================
 int i, j;
 
 if(opt_Ro<0)
  {
  log_Rmed=log_priorRo;
  }

  Lr=Par_bio(3);
  sr=Par_bio(4);

  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}


// genero la composicion de tallas del reclutamiento
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);

// genero una estructura inicial en torno a Z del primer año;
  Reclutas=mfexp(log_Rmed+log_desv_Rt);

// genero la poblacion en equilibrio virginal de LP;

  No=pre*exp(log_Rmed);
  for (int j=1;j<=5*nedades;j++){
  No=(No*exp(-1.*M))*T+pre*exp(log_Rmed); }
  
  SSBo    = sum(elem_prod(No*mfexp(-dt(1)*M),elem_prod(Wmed,msex)));
  alfa_sr = 4*h*exp(log_Rmed+0.5*square(sigmaR))/(5*h-1);//
  beta_sr = (1-h)*SSBo/(5*h-1);// Reclutamiento

// -----------------primer año
  Reclutas(1) = mfexp(log_Rmed+log_desv_Rt(1));
  Rpred(1)    = Reclutas(1);
// genero una estructura inicial en torno a Z del primer año;
  Neq=pre*Reclutas(1);
  for (j=1;j<=nedades;j++){
  Neq    = elem_prod(Neq,exp(-1.*Z(1)))*T+pre*exp(log_Rmed+log_desv_No(j));}
  N(1)   = Neq;
  NMD(1) = elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1))),msex);
  BD(1)  = sum(elem_prod(Wmed,NMD(1)));

// --------------------dinamica anual
  for (i=2;i<=ntime;i++){
  Reclutas(i) = mfexp(log_Rmed+log_desv_Rt(i));
  Rpred(i)    = Reclutas(i);
  
  if(i>minedad){
  Rpred(i)    = (alfa_sr*BD(i-minedad)/(beta_sr+BD(i-minedad)));
  Reclutas(i) = Rpred(i)*mfexp(log_desv_Rt(i)); }

  N(i)        = (elem_prod(N(i-1),S(i-1)))*T+pre*Reclutas(i);
  NMD(i)      = elem_prod(elem_prod(N(i),mfexp(-dt(1)*Z(i))),msex);
  BD(i)       = sum(elem_prod(Wmed,NMD(i)));
  } 

//========================================================================

FUNCTION Eval_biomasas

//========================================================================
 NMD=elem_prod(N,mfexp(-dt(1)*Z));
 NMD=elem_prod(NMD,outer_prod(Unos_year,msex));
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel);
 NVcru=elem_prod(elem_prod(N,mfexp(-dt(3)*(Z))),Selc);
 
 for(int i=1;i<=ntime;i++){
  NVcru(i)=elem_prod(elem_prod(N(i),mfexp(-dt_C(i)*(Z(i)))),Selc(i));
 }

// cout << "dt_C:" << dt_C << endl;exit(1);
// TEMP=outer_prod(ntallas,dt_C);
// NVcru=elem_prod(elem_prod(N,mfexp(outer_prod(dt_C,ntallas)*(-Z))),Selc);
// NVmph=elem_prod(elem_prod(N,mfexp(-dt(4)*(Z))),Selc);

// vectores de biomasas derivadas
 BD    = Wmed*trans(NMD);
 BMflo = Wmed*trans(NVflo);
 BMcru = Wmed*trans(NVcru);
 BMmph = Wmed*trans(NVmph);
 BT    = Wmed*trans(N);
 BDo   = sum(elem_prod(No,Wmed));
 RPRlp = BD/SSBo;
 RPR   = BD/(SSBo*0.55);
 Frpr  = elem_div(mfexp(log_F),Frms_bloque);
//========================================================================

FUNCTION Eval_capturas_predichas

//========================================================================

// matrices de capturas predichas por edad y año
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// vectores de desembarques predichos por año
 pred_Desemb=Wmed*trans(pred_Ctot);

// matrices de proporcion de capturas por talla y año
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));

 pobsc=elem_div(Ncru,outer_prod(rowsum(Ncru+1e-10),Unos_tallas));
 ppredc=elem_div(NVcru,outer_prod(rowsum(NVcru+1e-10),Unos_tallas));

 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);

 Lmed_predc=Tallas*trans(ppredc);
 Lmed_obsc=Tallas*trans(pobsc);

//========================================================================

FUNCTION Eval_indices

//========================================================================

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUE(i)=exp(log_qflo(j))*pow(BMflo(i),exp(log_b));}
       }
   }


   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesc;j++){
              if (yrs(i)>=yqbloquesc(j)){
                 pred_Bcru(i)=exp(log_qcru(j))*BMcru(i);}
       }
   }


  for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesmph;j++){
              if (yrs(i)>=yqbloquesmph(j)){
                 pred_mph(i)=exp(log_qmph(j))*BD(i);}
       }
   }


//========================================================================

FUNCTION Eval_logverosim

//========================================================================
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0; suma2=0; suma3=0; penalty=0;

 for (i=1;i<=ntime;i++)
 {
  if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}

  if (Bcru(i)>0){
    suma2+=square(log(Bcru(i)/pred_Bcru(i))*1/cv_index(3,i));}

   if (mph(i)>0){
    suma3+=square(log(mph(i)/pred_mph(i))*1/cv_index(4,i));}
 }

//========================================================================

FUNCTION Eval_funcion_objetivo

//========================================================================

 suma4=0; suma5=0; penalty=0;

 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*suma2;//Bcru
 likeval(3)=0.5*suma3;//mph

 likeval(4)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb

 for (int i=1;i<=ntime;i++){
 suma4+=-nm(i)*sum(elem_prod(pobs(i),log(ppred(i)+1e-10)));
 suma5+=-nmc(i)*sum(elem_prod(pobsc(i),log(ppredc(i)+1e-10)));
 }
 
 likeval(5)=suma4;//
 likeval(6)=suma5;//

// lognormal Ninicial y Reclutas
 if(active(log_desv_Rt)){
 likeval(7)=1./(2*square(sigmaR))*norm2(log_desv_Rt);}

 if(active(log_desv_No)){
 likeval(8)=1./(2*square(sigmaR))*norm2(log_desv_No);}

 if(active(log_Lr)){
 likeval(9)=1./(2*square(cv_Par_bio(3)))*square(log_Lr-log_Lr_prior);}

  //if (active(log_F)){
  //pF=1000*norm2(log_F-mean(log_F));}

  penalty(1)=0.5/square(cv_q(1))*norm2(log_qflo-log_prior_qf);
  penalty(2)=0.5/square(cv_q(2))*norm2(log_qcru-log_prior_qc);
  penalty(3)=0.5/square(cv_q(3))*norm2(log_qmph-log_prior_qmph);
  penalty(4)=0.5/square(0.4)*norm2(log_sigma1-log_s1prior);
  penalty(5)=0.5/square(0.4)*norm2(log_L50-log_L50prior);
  penalty(6)=0.5/square(cv_Par_bio(1))*square(log_Linf-log_Linf_prior);
  penalty(7)=0.5/square(cv_Par_bio(2))*square(log_k-log_k_prior);

  //f=festim*(sum(likeval)+sum(penalty)+pF);
  f=festim*(sum(likeval)+sum(penalty));
  
  if(last_phase){
  f=festim*(sum(likeval)+sum(penalty));}
 
//========================================================================

FUNCTION  Eval_CTP

//========================================================================

  for (int i=1;i<=npbr;i++){ // ciclo de PBR
  Npact=N(ntime);
  Np=N(ntime);
  Sp=S(ntime);

 // Fpbr=F(ntime)*pbr(i);// este usa multiplicador pbr
  Fpbr=Sel(ntime)*pbr(i);// este usa la selectividad x Fpbr
  Zpbr=Fpbr+M;
  NV      = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Npact));
  CTP     = elem_prod(NV,Wmed);
  Yact(i) = sum(CTP);

  for (int j=1;j<=ntime_sim;j++){ // ciclo de años

  if(j<=minedad){
 // Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*BD(ntime-minedad+1)/(beta_sr+BD(ntime-minedad+1)));} // Estima CTP con R_last
 if(opt_Ro<0)
  {
  log_Rmed=log_priorRo;
  }
  
  if(oprec==1){  
  Np=(elem_prod(Np,Sp))*T+pre*(mfexp(log_Rmed));} // Estima CTP con R_med
  if(oprec==2){  
  Np=(elem_prod(Np,Sp))*T+pre*Reclutas(11);} // Estima CTP con R_alto (2012)
  if(oprec==3){  
  Np=(elem_prod(Np,Sp))*T+pre*Reclutas(17);} // Estima CTP con R_bajo (2018)
  }
 if(j>minedad){
  Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*Bp(i,j-minedad)/(beta_sr+Bp(i,j-minedad)));} //

  Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dt(1)*Zpbr)),elem_prod(msex,Wmed)));
  
  NVp  = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  CTPp = elem_prod(NVp,Wmed);
  Yp(i,j)=sum(CTPp);
  Sp=exp(-1.*Zpbr);
  }}

//========================================================================

REPORT_SECTION

//========================================================================
 report << "Years" << endl;
 report << yrs << endl;
 report << "Bcru_obs" << endl;
 report << Bcru << endl;
 report << "Bcru_pred" << endl;
 report << pred_Bcru << endl;
 report << "CPUE_obs" << endl;
 report << CPUE << endl;
 report << "CPUE_pred" << endl;
 report << pred_CPUE << endl;
 report << "MPH_obs" << endl;
 report << mph << endl;
 report << "MPH_pred" << endl;
 report << pred_mph << endl;
 report << "Desemb_obs" << endl;
 report << Desemb << endl;
 report << "Desemb_pred" << endl;
 report << pred_Desemb << endl;
 report << "Lmf_obs" << endl;
 report << Lmed_obs << endl;
 report << "Lmf_pred" << endl;
 report << Lmed_pred << endl;
 report << "Lmc_obs" << endl;
 report << Lmed_obsc << endl;
 report << "Lmc_pred" << endl;
 report << Lmed_predc << endl;
 report << "Biomasa_desovante" << endl;
 report << BD << endl;
 report << "Biomasa_total" << endl;
 report << BT << endl;
 report << "Biomasa_explotable" << endl;
 report << BMflo << endl;
 report << "Reclutamiento" << endl;
 report << Reclutas<< endl;
 report << "Rpred" << endl;
 report << Rpred<< endl;
 report << "F" << endl;
 report << exp(log_F) << endl;
 report<<"Tallas"<<endl;
 report<<Tallas<<endl;
 report<<"Abundancia_talla"<<endl;
 report<<N<<endl;
 report<<"Selflo_talla"<<endl;
 report<<Sel<<endl;
 report<<"Selcru_talla"<<endl;
 report<<Selc<<endl;
 report << "Propfl_obs" << endl;
 report << pobs<< endl;
 report << "Propfl_pred" << endl;
 report << ppred<< endl;
 report << "Propcru_obs" << endl;
 report << pobsc<< endl;
 report << "Propcru_pred" << endl;
 report << ppredc<< endl;
 report << "BD_virgen_anual" << endl;
 report << BDo << endl;
 report << "BD_virgen_LP" << endl;
 report << SSBo << endl;
 report << "Reduccion_LP " << endl;
 report << RPRlp << endl;
 report << "Talla_media_por_grupo" << endl;
 report << Lesp << endl;
 report <<  "desvest_por_grupo" << endl;
 report << sigmaL << endl;
 report << "Fun_rec_talla" << endl;
 report << pre<< endl;
 report << "MatrizTrans" << endl;
 report << T << endl;
 report << "bCPUE  Lr  Sr  beta  h " << endl;
 report << exp(log_b)<<" "<<exp(log_Lr)<<" "<<exp(log_sr)<<" "<<exp(log_beta)<<" "<<h<< endl;
 report << "pred_Ctot" << endl;
 report << pred_Ctot << endl;

//-------------------------------------------------------------------
// ESTIMA nm y CV

  suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppred(i),1-ppred(i)));
      suma2=norm2(pobs(i)-ppred(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report << "nm_flota_cru" <<endl;
 report <<pow(nm1,1/cuenta1)<< endl;


 suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppredc(i),1-ppredc(i)));
      suma2=norm2(pobsc(i)-ppredc(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report <<pow(nm1,1/cuenta1)<< endl;
 report << "BD_proy" << endl; //biomasa desovante proyectada para cada multiplo de Flast" 
 report << Bp << endl;
 report << "Capt_proy" << endl; // Capturas proyectadas para cada Fpbr
 report << Yp << endl;
 report << "Captura_act" << endl;  //Captura proyectadas año en curso
 report << Yact << endl;
 report<<"Likeval CPUE_BCru_Bmph_Desemb_pfFlota_pfCru_dev_R_devNo_LR"<<endl;
 report << likeval << endl;
 report << "Prioris" << endl;
 report << penalty << endl;
 report << "Rec_proy" << endl; //Reclutamientos proyectadas para cada Fpbr
 report << Rpp << endl;
 report << "N_actual" << endl; // abundancia del ultimo año
 report << Npact << endl;
 report << "N_proyect" << endl; // abundancia del ultimo año
 report << Np << endl;
 report << "Rproy1" << endl; //Rec
 report << mfexp(log_Rmed) << endl;
 report << "Rproy2" << endl; //Rec
 report << Reclutas(11)<< endl;
 report << "Rproy3" << endl; //Rec
 report << Reclutas(17)<< endl;
 
 report << "NV" << endl;
 report << NV << endl;
 report << "CTP" << endl;
 report << CTP << endl;
 report << "NVp" << endl;
 report << NVp << endl;
 report << "CTPp" << endl;
 report << CTPp << endl;

 
FINAL_SECTION

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;


