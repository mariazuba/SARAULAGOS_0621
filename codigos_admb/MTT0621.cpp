#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <MTT0621.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  ntime.allocate("ntime");
  nedades.allocate("nedades");
  ntallas.allocate("ntallas");
  mdatos.allocate(1,ntime,1,13,"mdatos");
  Tallas.allocate(1,ntallas,"Tallas");
  Ctot.allocate(1,ntime,1,ntallas,"Ctot");
  Ncru.allocate(1,ntime,1,ntallas,"Ncru");
  msex.allocate(1,ntallas,"msex");
  Wmed.allocate(1,ntallas,"Wmed");
  sigmaR.allocate("sigmaR");
  dt.allocate(1,3,"dt");
  Par_bio.allocate(1,7,"Par_bio");
  cv_Par_bio.allocate(1,7,"cv_Par_bio");
  minedad.allocate("minedad");
  bprior.allocate("bprior");
 log_Linf_prior = log(Par_bio(1));
 log_k_prior = log(Par_bio(2));
 log_Lr_prior = log(Par_bio(3));
 log_sr_prior = log(Par_bio(4));
 log_beta_prior= log(Par_bio(5));
 log_b_prior = log(bprior);
  L50prior.allocate("L50prior");
  s1prior.allocate("s1prior");
  s2prior.allocate("s2prior");
  opt_sel2.allocate("opt_sel2");
  opt_sel3.allocate("opt_sel3");
 log_L50prior = log(L50prior);
 log_s1prior = log(s1prior);
 log_s2prior = log(s2prior);
  nbloques1.allocate("nbloques1");
  ybloques1.allocate(1,nbloques1,"ybloques1");
  nbloques2.allocate("nbloques2");
  ybloques2.allocate(1,nbloques2,"ybloques2");
  nqbloques.allocate("nqbloques");
  yqbloques.allocate(1,nqbloques,"yqbloques");
  nqbloquesc.allocate("nqbloquesc");
  yqbloquesc.allocate(1,nqbloquesc,"yqbloquesc");
  nqbloquesmph.allocate("nqbloquesmph");
  yqbloquesmph.allocate(1,nqbloquesmph,"yqbloquesmph");
  prior_qf.allocate("prior_qf");
  prior_qc.allocate("prior_qc");
  prior_qmph.allocate("prior_qmph");
log_prior_qf=log(prior_qf);
log_prior_qc=log(prior_qc); 
log_prior_qmph =log(prior_qmph);
  cv_q.allocate(1,3,"cv_q");
  opt_qf.allocate("opt_qf");
  opt_qc.allocate("opt_qc");
  opt_qmph.allocate("opt_qmph");
  opt1_fase.allocate("opt1_fase");
  opt2_fase.allocate("opt2_fase");
  opt_Lr.allocate("opt_Lr");
  opt_sr.allocate("opt_sr");
  opt_beta.allocate("opt_beta");
  opt_vb1.allocate("opt_vb1");
  opt_vb2.allocate("opt_vb2");
  opt_Ro.allocate("opt_Ro");
  log_priorRo.allocate("log_priorRo");
  opt_F.allocate("opt_F");
  opt_devRt.allocate("opt_devRt");
  opt_devNo.allocate("opt_devNo");
  opt_bpow.allocate("opt_bpow");
  npbr.allocate("npbr");
  pbr.allocate(1,npbr,"pbr");
  ntime_sim.allocate("ntime_sim");
  festim.allocate("festim");
  oprec.allocate("oprec");
}

void model_parameters::initializationfunction(void)
{
  log_Rmed.set_initial_value(8.8);
  log_Lr.set_initial_value(log_Lr_prior);
  log_sr.set_initial_value(log_sr_prior);
  log_L50.set_initial_value(log_L50prior);
  log_L50c.set_initial_value(log_L50prior);
  log_sigma1.set_initial_value(log_s1prior);
  log_sigma2.set_initial_value(log_s2prior);
  log_sigma1c.set_initial_value(log_s1prior);
  log_sigma2c.set_initial_value(log_s2prior);
  log_b.set_initial_value(log_b_prior);
  log_beta.set_initial_value(log_beta_prior);
  log_k.set_initial_value(log_k_prior);
  log_Linf.set_initial_value(log_Linf_prior);
  log_qflo.set_initial_value(log_prior_qf);
  log_qcru.set_initial_value(log_prior_qc);
  log_qmph.set_initial_value(log_prior_qmph);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_L50.allocate(1,nbloques1,opt1_fase,"log_L50");
  log_sigma1.allocate(1,nbloques1,opt1_fase,"log_sigma1");
  log_sigma2.allocate(1,nbloques1,opt_sel2,"log_sigma2");
  log_L50c.allocate(1,nbloques2,opt2_fase,"log_L50c");
  log_sigma1c.allocate(1,nbloques2,opt2_fase,"log_sigma1c");
  log_sigma2c.allocate(1,nbloques2,opt_sel3,"log_sigma2c");
  log_Rmed.allocate(opt_Ro,"log_Rmed");
  log_desv_Rt.allocate(1,ntime,-10,10,opt_devRt,"log_desv_Rt");
  log_desv_No.allocate(1,nedades,-10,10,opt_devNo,"log_desv_No");
  log_F.allocate(1,ntime,-20,0.7,opt_F,"log_F");
  log_qflo.allocate(1,nqbloques,opt_qf,"log_qflo");
  log_qcru.allocate(1,nqbloquesc,opt_qc,"log_qcru");
  log_qmph.allocate(1,nqbloquesmph,opt_qmph,"log_qmph");
  log_b.allocate(opt_bpow,"log_b");
  log_Lr.allocate(opt_Lr,"log_Lr");
  log_sr.allocate(opt_sr,"log_sr");
  log_beta.allocate(opt_beta,"log_beta");
  log_Linf.allocate(opt_vb1,"log_Linf");
  log_k.allocate(opt_vb2,"log_k");
  BMflo.allocate(1,ntime,"BMflo");
  #ifndef NO_AD_INITIALIZE
    BMflo.initialize();
  #endif
  BMcru.allocate(1,ntime,"BMcru");
  #ifndef NO_AD_INITIALIZE
    BMcru.initialize();
  #endif
  Brec.allocate(1,ntime,"Brec");
  #ifndef NO_AD_INITIALIZE
    Brec.initialize();
  #endif
  BMmph.allocate(1,ntime,"BMmph");
  #ifndef NO_AD_INITIALIZE
    BMmph.initialize();
  #endif
  pred_CPUE.allocate(1,ntime,"pred_CPUE");
  pred_Bcru.allocate(1,ntime,"pred_Bcru");
  pred_Desemb.allocate(1,ntime,"pred_Desemb");
  likeval.allocate(1,9,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  Neq.allocate(1,ntallas,"Neq");
  #ifndef NO_AD_INITIALIZE
    Neq.initialize();
  #endif
  Rpred.allocate(1,ntime,"Rpred");
  Unos_edad.allocate(1,nedades,"Unos_edad");
  #ifndef NO_AD_INITIALIZE
    Unos_edad.initialize();
  #endif
  Unos_year.allocate(1,ntime,"Unos_year");
  #ifndef NO_AD_INITIALIZE
    Unos_year.initialize();
  #endif
  Unos_tallas.allocate(1,ntallas,"Unos_tallas");
  #ifndef NO_AD_INITIALIZE
    Unos_tallas.initialize();
  #endif
  delta.allocate(1,ntallas,"delta");
  #ifndef NO_AD_INITIALIZE
    delta.initialize();
  #endif
  Lesp.allocate(1,ntallas,"Lesp");
  #ifndef NO_AD_INITIALIZE
    Lesp.initialize();
  #endif
  sigmaL.allocate(1,ntallas,"sigmaL");
  #ifndef NO_AD_INITIALIZE
    sigmaL.initialize();
  #endif
  pre.allocate(1,ntallas,"pre");
  #ifndef NO_AD_INITIALIZE
    pre.initialize();
  #endif
  mu_edad.allocate(1,nedades,"mu_edad");
  #ifndef NO_AD_INITIALIZE
    mu_edad.initialize();
  #endif
  sigma_edad.allocate(1,nedades,"sigma_edad");
  #ifndef NO_AD_INITIALIZE
    sigma_edad.initialize();
  #endif
  BDo.allocate(1,ntime,"BDo");
  #ifndef NO_AD_INITIALIZE
    BDo.initialize();
  #endif
  No.allocate(1,ntallas,"No");
  #ifndef NO_AD_INITIALIZE
    No.initialize();
  #endif
  prior.allocate(1,7,"prior");
  #ifndef NO_AD_INITIALIZE
    prior.initialize();
  #endif
  yrs.allocate(1,ntime,"yrs");
  #ifndef NO_AD_INITIALIZE
    yrs.initialize();
  #endif
  Desemb.allocate(1,ntime,"Desemb");
  #ifndef NO_AD_INITIALIZE
    Desemb.initialize();
  #endif
  CPUE.allocate(1,ntime,"CPUE");
  #ifndef NO_AD_INITIALIZE
    CPUE.initialize();
  #endif
  Bcru.allocate(1,ntime,"Bcru");
  #ifndef NO_AD_INITIALIZE
    Bcru.initialize();
  #endif
  dt_C.allocate(1,ntime,"dt_C");
  #ifndef NO_AD_INITIALIZE
    dt_C.initialize();
  #endif
  Frms_bloque.allocate(1,ntime,"Frms_bloque");
  #ifndef NO_AD_INITIALIZE
    Frms_bloque.initialize();
  #endif
  mph.allocate(1,ntime,"mph");
  #ifndef NO_AD_INITIALIZE
    mph.initialize();
  #endif
  Lmed_obs.allocate(1,ntime,"Lmed_obs");
  #ifndef NO_AD_INITIALIZE
    Lmed_obs.initialize();
  #endif
  Lmed_pred.allocate(1,ntime,"Lmed_pred");
  Lmed_obsc.allocate(1,ntime,"Lmed_obsc");
  #ifndef NO_AD_INITIALIZE
    Lmed_obsc.initialize();
  #endif
  Lmed_predc.allocate(1,ntime,"Lmed_predc");
  edades.allocate(1,nedades,"edades");
  #ifndef NO_AD_INITIALIZE
    edades.initialize();
  #endif
  Reclutas.allocate(1,ntime,"Reclutas");
  nm.allocate(1,ntime,"nm");
  #ifndef NO_AD_INITIALIZE
    nm.initialize();
  #endif
  nmc.allocate(1,ntime,"nmc");
  #ifndef NO_AD_INITIALIZE
    nmc.initialize();
  #endif
  penalty.allocate(1,7,"penalty");
  #ifndef NO_AD_INITIALIZE
    penalty.initialize();
  #endif
  cv_index.allocate(1,4,1,ntime,"cv_index");
  #ifndef NO_AD_INITIALIZE
    cv_index.initialize();
  #endif
  S1.allocate(1,nbloques1,1,ntallas,"S1");
  #ifndef NO_AD_INITIALIZE
    S1.initialize();
  #endif
  S2.allocate(1,nbloques2,1,ntallas,"S2");
  #ifndef NO_AD_INITIALIZE
    S2.initialize();
  #endif
  Sel.allocate(1,ntime,1,ntallas,"Sel");
  #ifndef NO_AD_INITIALIZE
    Sel.initialize();
  #endif
  Selc.allocate(1,ntime,1,ntallas,"Selc");
  #ifndef NO_AD_INITIALIZE
    Selc.initialize();
  #endif
  F.allocate(1,ntime,1,ntallas,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,ntime,1,ntallas,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,ntime,1,ntallas,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,ntime,1,ntallas,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  NM.allocate(1,ntime,1,ntallas,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  NMD.allocate(1,ntime,1,ntallas,"NMD");
  #ifndef NO_AD_INITIALIZE
    NMD.initialize();
  #endif
  NDv.allocate(1,ntime,1,ntallas,"NDv");
  #ifndef NO_AD_INITIALIZE
    NDv.initialize();
  #endif
  Nrec.allocate(1,ntime,1,ntallas,"Nrec");
  #ifndef NO_AD_INITIALIZE
    Nrec.initialize();
  #endif
  NVflo.allocate(1,ntime,1,ntallas,"NVflo");
  #ifndef NO_AD_INITIALIZE
    NVflo.initialize();
  #endif
  NVcru.allocate(1,ntime,1,ntallas,"NVcru");
  #ifndef NO_AD_INITIALIZE
    NVcru.initialize();
  #endif
  NVmph.allocate(1,ntime,1,ntallas,"NVmph");
  #ifndef NO_AD_INITIALIZE
    NVmph.initialize();
  #endif
  TEMP.allocate(1,ntime,1,ntallas,"TEMP");
  #ifndef NO_AD_INITIALIZE
    TEMP.initialize();
  #endif
  pred_Ctot.allocate(1,ntime,1,ntallas,"pred_Ctot");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot.initialize();
  #endif
  pobs.allocate(1,ntime,1,ntallas,"pobs");
  #ifndef NO_AD_INITIALIZE
    pobs.initialize();
  #endif
  ppred.allocate(1,ntime,1,ntallas,"ppred");
  #ifndef NO_AD_INITIALIZE
    ppred.initialize();
  #endif
  pobsc.allocate(1,ntime,1,ntallas,"pobsc");
  #ifndef NO_AD_INITIALIZE
    pobsc.initialize();
  #endif
  ppredc.allocate(1,ntime,1,ntallas,"ppredc");
  #ifndef NO_AD_INITIALIZE
    ppredc.initialize();
  #endif
  T.allocate(1,ntallas,1,ntallas,"T");
  #ifndef NO_AD_INITIALIZE
    T.initialize();
  #endif
  Nv.allocate(1,ntime,1,nedades,"Nv");
  #ifndef NO_AD_INITIALIZE
    Nv.initialize();
  #endif
  NMDv.allocate(1,ntime,1,ntallas,"NMDv");
  #ifndef NO_AD_INITIALIZE
    NMDv.initialize();
  #endif
  suma1.allocate("suma1");
  #ifndef NO_AD_INITIALIZE
  suma1.initialize();
  #endif
  suma2.allocate("suma2");
  #ifndef NO_AD_INITIALIZE
  suma2.initialize();
  #endif
  suma3.allocate("suma3");
  #ifndef NO_AD_INITIALIZE
  suma3.initialize();
  #endif
  suma4.allocate("suma4");
  #ifndef NO_AD_INITIALIZE
  suma4.initialize();
  #endif
  suma5.allocate("suma5");
  #ifndef NO_AD_INITIALIZE
  suma5.initialize();
  #endif
  So.allocate("So");
  #ifndef NO_AD_INITIALIZE
  So.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  Linfh.allocate("Linfh");
  #ifndef NO_AD_INITIALIZE
  Linfh.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  Lr.allocate("Lr");
  #ifndef NO_AD_INITIALIZE
  Lr.initialize();
  #endif
  sr.allocate("sr");
  #ifndef NO_AD_INITIALIZE
  sr.initialize();
  #endif
  Lm.allocate("Lm");
  #ifndef NO_AD_INITIALIZE
  Lm.initialize();
  #endif
  Rm.allocate("Rm");
  #ifndef NO_AD_INITIALIZE
  Rm.initialize();
  #endif
  h.allocate("h");
  #ifndef NO_AD_INITIALIZE
  h.initialize();
  #endif
  BDp.allocate("BDp");
  #ifndef NO_AD_INITIALIZE
  BDp.initialize();
  #endif
  Npplus.allocate("Npplus");
  #ifndef NO_AD_INITIALIZE
  Npplus.initialize();
  #endif
  Bp_anch.allocate("Bp_anch");
  #ifndef NO_AD_INITIALIZE
  Bp_anch.initialize();
  #endif
  nm1.allocate("nm1");
  #ifndef NO_AD_INITIALIZE
  nm1.initialize();
  #endif
  cuenta1.allocate("cuenta1");
  #ifndef NO_AD_INITIALIZE
  cuenta1.initialize();
  #endif
  alfa_sr.allocate("alfa_sr");
  #ifndef NO_AD_INITIALIZE
  alfa_sr.initialize();
  #endif
  beta_sr.allocate("beta_sr");
  #ifndef NO_AD_INITIALIZE
  beta_sr.initialize();
  #endif
  pF.allocate("pF");
  #ifndef NO_AD_INITIALIZE
  pF.initialize();
  #endif
  Npact.allocate(1,ntallas,"Npact");
  #ifndef NO_AD_INITIALIZE
    Npact.initialize();
  #endif
  Np.allocate(1,ntallas,"Np");
  #ifndef NO_AD_INITIALIZE
    Np.initialize();
  #endif
  Zpbr.allocate(1,ntallas,"Zpbr");
  #ifndef NO_AD_INITIALIZE
    Zpbr.initialize();
  #endif
  Fpbr.allocate(1,ntallas,"Fpbr");
  Sp.allocate(1,ntallas,"Sp");
  #ifndef NO_AD_INITIALIZE
    Sp.initialize();
  #endif
  Bp.allocate(1,npbr,1,ntime_sim,"Bp");
  NV.allocate(1,ntallas,"NV");
  #ifndef NO_AD_INITIALIZE
    NV.initialize();
  #endif
  CTP.allocate(1,ntallas,"CTP");
  #ifndef NO_AD_INITIALIZE
    CTP.initialize();
  #endif
  Yact.allocate(1,npbr,"Yact");
  NVp.allocate(1,ntallas,"NVp");
  #ifndef NO_AD_INITIALIZE
    NVp.initialize();
  #endif
  CTPp.allocate(1,ntallas,"CTPp");
  #ifndef NO_AD_INITIALIZE
    CTPp.initialize();
  #endif
  Yp.allocate(1,npbr,1,ntime_sim,"Yp");
  Rpp.allocate(1,npbr,1,ntime_sim,"Rpp");
  #ifndef NO_AD_INITIALIZE
    Rpp.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  BD.allocate(1,ntime,"BD");
  BT.allocate(1,ntime,"BT");
  RPRlp.allocate(1,ntime,"RPRlp");
  RPR.allocate(1,ntime,"RPR");
  Frpr.allocate(1,ntime,"Frpr");
  pred_mph.allocate(1,ntime,"pred_mph");
  SSBo.allocate("SSBo");
  SPRFo.allocate("SPRFo");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
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
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e-1,1.e-01,1.e-03,1e-3,1e-5}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{100,100,200,3000,3500}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::userfunction(void)
{
  f =0.0;
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
}

void model_parameters::Eval_Trans_talla_talla(void)
{
  Linf=exp(log_Linf);
  k=exp(log_k);
  beta=exp(log_beta);
  if(active(log_beta)){beta=mfexp(log_beta);}
 int i, j;
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
}

void model_parameters::Eval_selectividad(void)
{
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
}

void model_parameters::Eval_mortalidades(void)
{
 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_tallas));
 Z=F+M;
 S=mfexp(-1.0*Z);
}

void model_parameters::Eval_abundancia(void)
{
 int i, j;
 if(opt_Ro<0)
  {
  log_Rmed=log_priorRo;
  }
  Lr=Par_bio(3);
  sr=Par_bio(4);
  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);
  Reclutas=mfexp(log_Rmed+log_desv_Rt);
  No=pre*exp(log_Rmed);
  for (int j=1;j<=5*nedades;j++){
  No=(No*exp(-1.*M))*T+pre*exp(log_Rmed); }
  SSBo    = sum(elem_prod(No*mfexp(-dt(1)*M),elem_prod(Wmed,msex)));
  SPRFo    = sum(elem_prod(No*mfexp(-dt(1)*M),msex));//para RAM_legacy
  alfa_sr = 4*h*exp(log_Rmed+0.5*square(sigmaR))/(5*h-1);//
  beta_sr = (1-h)*SSBo/(5*h-1);// Reclutamiento
  Reclutas(1) = mfexp(log_Rmed+log_desv_Rt(1));
  Rpred(1)    = Reclutas(1);
  Neq=pre*Reclutas(1);
  for (j=1;j<=nedades;j++){
  Neq    = elem_prod(Neq,exp(-1.*Z(1)))*T+pre*exp(log_Rmed+log_desv_No(j));}
  N(1)   = Neq;
  NMD(1) = elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1))),msex);
  BD(1)  = sum(elem_prod(Wmed,NMD(1)));
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
}

void model_parameters::Eval_biomasas(void)
{
 NMD=elem_prod(N,mfexp(-dt(1)*Z));
 NMD=elem_prod(NMD,outer_prod(Unos_year,msex));
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel);
 NVcru=elem_prod(elem_prod(N,mfexp(-dt(3)*(Z))),Selc);
 for(int i=1;i<=ntime;i++){
  NVcru(i)=elem_prod(elem_prod(N(i),mfexp(-dt_C(i)*(Z(i)))),Selc(i));
 }
 BD    = Wmed*trans(NMD);
 BMflo = Wmed*trans(NVflo);
 BMcru = Wmed*trans(NVcru);
 BMmph = Wmed*trans(NVmph);
 BT    = Wmed*trans(N);
 BDo   = sum(elem_prod(No,Wmed));
 RPRlp = BD/SSBo;
 RPR   = BD/(SSBo*0.55);
 Frpr  = elem_div(mfexp(log_F),Frms_bloque);
}

void model_parameters::Eval_capturas_predichas(void)
{
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
 pred_Desemb=Wmed*trans(pred_Ctot);
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));
 pobsc=elem_div(Ncru,outer_prod(rowsum(Ncru+1e-10),Unos_tallas));
 ppredc=elem_div(NVcru,outer_prod(rowsum(NVcru+1e-10),Unos_tallas));
 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);
 Lmed_predc=Tallas*trans(ppredc);
 Lmed_obsc=Tallas*trans(pobsc);
}

void model_parameters::Eval_indices(void)
{
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
}

void model_parameters::Eval_logverosim(void)
{
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
}

void model_parameters::Eval_funcion_objetivo(void)
{
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
}

void model_parameters::Eval_CTP(void)
{
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
 report<<"pred_Ctot"<<endl;
 report<<pred_Ctot<<endl;
 report<<"NVcru"<<endl;
 report<<NVcru<<endl;
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
 report << "NBD_virgen_LP" << endl;
 report <<  SPRFo << endl;
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
}

void model_parameters::final_calcs()
{
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
