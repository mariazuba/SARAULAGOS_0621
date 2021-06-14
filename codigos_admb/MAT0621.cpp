#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 #include <limits>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 double eps=std::numeric_limits<double>::epsilon();
 #include <string.h>
 #undef depur
 #undef depuro
 #define depur(object) cout << #object "\n" << object << endl; 
 #define depuro(object) cout << #object "\n" << object << endl; exit(1);
 //#########################################################
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
#include <MAT0621.htp>

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
  edades.allocate(1,nedades,"edades");
  Tallas.allocate(1,ntallas,"Tallas");
  data.allocate(1,ntime,1,13,"data");
  Ctot.allocate(1,ntime,1,ntallas,"Ctot");
  Ncru.allocate(1,ntime,1,ntallas,"Ncru");
  msex.allocate(1,ntallas,"msex");
  Wmed.allocate(1,ntallas,"Wmed");
  cvar.allocate(1,2,"cvar");
  dt.allocate(1,3,"dt");
  Par_bio.allocate(1,7,"Par_bio");
  cv_priors.allocate(1,7,"cv_priors");
  fases_bio.allocate(1,7,"fases_bio");
 log_Linfprior   = log(Par_bio(1));
 log_kprior      = log(Par_bio(2));
 log_Loprior     = log(Par_bio(3));
 log_alfa_prior  = log(Par_bio(4));
 log_beta_prior  = log(Par_bio(5));
 log_M_prior     = log(Par_bio(6));
 log_h_prior     = log(Par_bio(7));
 opt_Linf  = fases_bio(1);
 opt_k     = fases_bio(2);
 opt_Lo    = fases_bio(3);
 opt_alfa  = fases_bio(4);
 opt_beta  = fases_bio(5);
 opt_M     = fases_bio(6);
 opt_h     = fases_bio(7);
  prior_qf.allocate("prior_qf");
  prior_qc.allocate("prior_qc");
  prior_qmph.allocate("prior_qmph");
  cv_qf.allocate("cv_qf");
  cv_qcru.allocate("cv_qcru");
  cv_qmph.allocate("cv_qmph");
 log_prior_qf = log(prior_qf);
 log_prior_qc = log(prior_qc);
 log_prior_qmph = log(prior_qmph);
  parf.allocate(1,2,"parf");
  parc.allocate(1,2,"parc");
 log_L50priorf = log(parf(1));
 log_spriorf   = log(parf(2));
 log_L50priorc = log(parc(1));
 log_spriorc   = log(parc(2));
  nbloques1.allocate("nbloques1");
  ybloques1.allocate(1,nbloques1,"ybloques1");
  nbloques2.allocate("nbloques2");
  ybloques2.allocate(1,nbloques2,"ybloques2");
  nqbloques.allocate("nqbloques");
  yqbloques.allocate(1,nqbloques,"yqbloques");
  nqbloquesc.allocate("nqbloquesc");
  yqbloquesc.allocate(1,nqbloquesc,"yqbloquesc");
  opt_qf.allocate("opt_qf");
  opt_qc.allocate("opt_qc");
  opt_qmph.allocate("opt_qmph");
  optSf_fase.allocate("optSf_fase");
  optS1_fase.allocate("optS1_fase");
  opt_F.allocate("opt_F");
  opt_devRt.allocate("opt_devRt");
  opt_devNo.allocate("opt_devNo");
  opt_Reclu.allocate("opt_Reclu");
  opt_Ro.allocate("opt_Ro");
  log_priorRo.allocate("log_priorRo");
  opt_Fpbr.allocate("opt_Fpbr");
  npbr.allocate("npbr");
  pbr.allocate(1,npbr,"pbr");
  ntime_sim.allocate("ntime_sim");
  pR.allocate("pR");
  oprec.allocate("oprec");
  opt_sim.allocate("opt_sim");
  bprior.allocate("bprior");
 log_b_prior = log(bprior);
  opt_bpow.allocate("opt_bpow");
}

void model_parameters::initializationfunction(void)
{
  log_Ro.set_initial_value(8.5);
  log_Linf.set_initial_value(log_Linfprior);
  log_k.set_initial_value(log_kprior);
  log_Lo.set_initial_value(log_Loprior);
  log_alfa.set_initial_value(log_alfa_prior);
  log_beta.set_initial_value(log_beta_prior);
  log_L50f.set_initial_value(log_L50priorf);
  log_rangof.set_initial_value(log_spriorf);
  log_L50c.set_initial_value(log_L50priorc);
  log_rangoc.set_initial_value(log_spriorc);
  log_M.set_initial_value(log_M_prior);
  log_qflo.set_initial_value(log_prior_qf);
  log_qcru.set_initial_value(log_prior_qc);
  log_b.set_initial_value(log_b_prior);
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
  log_L50f.allocate(1,nbloques1,optSf_fase,"log_L50f");
  log_rangof.allocate(1,nbloques1,optSf_fase,"log_rangof");
  log_L50c.allocate(1,nbloques2,optS1_fase,"log_L50c");
  log_rangoc.allocate(1,nbloques2,optS1_fase,"log_rangoc");
  log_Ro.allocate(opt_Ro,"log_Ro");
  dev_log_Ro.allocate(1,ntime,-10,10,opt_devRt,"dev_log_Ro");
  dev_log_No.allocate(1,nedades,-10,10,opt_devNo,"dev_log_No");
  log_F.allocate(1,ntime,-20,0.7,opt_F,"log_F");
  log_qflo.allocate(1,nqbloques,opt_qf,"log_qflo");
  log_qcru.allocate(1,nqbloquesc,opt_qc,"log_qcru");
  log_b.allocate(opt_bpow,"log_b");
  log_Linf.allocate(opt_Linf,"log_Linf");
  log_k.allocate(opt_k,"log_k");
  log_Lo.allocate(opt_Lo,"log_Lo");
  log_M.allocate(opt_M,"log_M");
  log_alfa.allocate(opt_alfa,"log_alfa");
  log_beta.allocate(opt_beta,"log_beta");
  log_Fref.allocate(1,npbr,opt_Fpbr,"log_Fref");
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
  mph.allocate(1,ntime,"mph");
  #ifndef NO_AD_INITIALIZE
    mph.initialize();
  #endif
  dt_C.allocate(1,ntime,"dt_C");
  #ifndef NO_AD_INITIALIZE
    dt_C.initialize();
  #endif
  Frms_bloque.allocate(1,ntime,"Frms_bloque");
  #ifndef NO_AD_INITIALIZE
    Frms_bloque.initialize();
  #endif
  cv_index.allocate(1,4,1,ntime,"cv_index");
  #ifndef NO_AD_INITIALIZE
    cv_index.initialize();
  #endif
  nm.allocate(1,2,1,ntime,"nm");
  #ifndef NO_AD_INITIALIZE
    nm.initialize();
  #endif
  Unos_edad.allocate(1,nedades,"Unos_edad");
  #ifndef NO_AD_INITIALIZE
    Unos_edad.initialize();
  #endif
  Unos_anos.allocate(1,ntime,"Unos_anos");
  #ifndef NO_AD_INITIALIZE
    Unos_anos.initialize();
  #endif
  Unos_tallas.allocate(1,ntallas,"Unos_tallas");
  #ifndef NO_AD_INITIALIZE
    Unos_tallas.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  Lo.allocate("Lo");
  #ifndef NO_AD_INITIALIZE
  Lo.initialize();
  #endif
  mu_edad.allocate(1,nedades,"mu_edad");
  #ifndef NO_AD_INITIALIZE
    mu_edad.initialize();
  #endif
  sigma_edad.allocate(1,nedades,"sigma_edad");
  #ifndef NO_AD_INITIALIZE
    sigma_edad.initialize();
  #endif
  Prob_talla.allocate(1,nedades,1,ntallas,"Prob_talla");
  #ifndef NO_AD_INITIALIZE
    Prob_talla.initialize();
  #endif
  S1.allocate(1,nbloques1,1,nedades,"S1");
  #ifndef NO_AD_INITIALIZE
    S1.initialize();
  #endif
  Sel_f.allocate(1,ntime,1,nedades,"Sel_f");
  #ifndef NO_AD_INITIALIZE
    Sel_f.initialize();
  #endif
  S2.allocate(1,nbloques2,1,nedades,"S2");
  #ifndef NO_AD_INITIALIZE
    S2.initialize();
  #endif
  Sel_c.allocate(1,ntime,1,nedades,"Sel_c");
  #ifndef NO_AD_INITIALIZE
    Sel_c.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  F.allocate(1,ntime,1,nedades,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,ntime,1,nedades,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,ntime,1,nedades,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  No.allocate(1,nedades,"No");
  #ifndef NO_AD_INITIALIZE
    No.initialize();
  #endif
  SSBo.allocate("SSBo");
  phi.allocate("phi");
  #ifndef NO_AD_INITIALIZE
  phi.initialize();
  #endif
  h.allocate("h");
  #ifndef NO_AD_INITIALIZE
  h.initialize();
  #endif
  alfa.allocate("alfa");
  #ifndef NO_AD_INITIALIZE
  alfa.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  Neq.allocate(1,nedades,"Neq");
  #ifndef NO_AD_INITIALIZE
    Neq.initialize();
  #endif
  N.allocate(1,ntime,1,nedades,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  Ntallas.allocate(1,ntime,1,ntallas,"Ntallas");
  #ifndef NO_AD_INITIALIZE
    Ntallas.initialize();
  #endif
  BD.allocate(1,ntime,"BD");
  Rpred.allocate(1,ntime,"Rpred");
  #ifndef NO_AD_INITIALIZE
    Rpred.initialize();
  #endif
  Reclutas.allocate(1,ntime,"Reclutas");
  Lf_obs.allocate(1,ntime,"Lf_obs");
  #ifndef NO_AD_INITIALIZE
    Lf_obs.initialize();
  #endif
  Lf_pred.allocate(1,ntime,"Lf_pred");
  Lc_obs.allocate(1,ntime,"Lc_obs");
  #ifndef NO_AD_INITIALIZE
    Lc_obs.initialize();
  #endif
  Lc_pred.allocate(1,ntime,"Lc_pred");
  Nv.allocate(1,ntime,1,nedades,"Nv");
  #ifndef NO_AD_INITIALIZE
    Nv.initialize();
  #endif
  NDv.allocate(1,ntime,1,ntallas,"NDv");
  #ifndef NO_AD_INITIALIZE
    NDv.initialize();
  #endif
  BDo.allocate(1,ntime,"BDo");
  #ifndef NO_AD_INITIALIZE
    BDo.initialize();
  #endif
  RPR.allocate(1,ntime,"RPR");
  RPRlp.allocate(1,ntime,"RPRlp");
  NMD.allocate(1,ntime,1,ntallas,"NMD");
  #ifndef NO_AD_INITIALIZE
    NMD.initialize();
  #endif
  NVflo.allocate(1,ntime,1,ntallas,"NVflo");
  #ifndef NO_AD_INITIALIZE
    NVflo.initialize();
  #endif
  NVcru.allocate(1,ntime,1,ntallas,"NVcru");
  #ifndef NO_AD_INITIALIZE
    NVcru.initialize();
  #endif
  BMflo.allocate(1,ntime,"BMflo");
  #ifndef NO_AD_INITIALIZE
    BMflo.initialize();
  #endif
  BMcru.allocate(1,ntime,"BMcru");
  #ifndef NO_AD_INITIALIZE
    BMcru.initialize();
  #endif
  BT.allocate(1,ntime,"BT");
  BV.allocate(1,ntime,"BV");
  pred_Ctot_a.allocate(1,ntime,1,nedades,"pred_Ctot_a");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot_a.initialize();
  #endif
  pred_Ctot.allocate(1,ntime,1,ntallas,"pred_Ctot");
  #ifndef NO_AD_INITIALIZE
    pred_Ctot.initialize();
  #endif
  pred_Desemb.allocate(1,ntime,"pred_Desemb");
  #ifndef NO_AD_INITIALIZE
    pred_Desemb.initialize();
  #endif
  pobs_f.allocate(1,ntime,1,ntallas,"pobs_f");
  #ifndef NO_AD_INITIALIZE
    pobs_f.initialize();
  #endif
  ppred_f.allocate(1,ntime,1,ntallas,"ppred_f");
  #ifndef NO_AD_INITIALIZE
    ppred_f.initialize();
  #endif
  pobsc.allocate(1,ntime,1,ntallas,"pobsc");
  #ifndef NO_AD_INITIALIZE
    pobsc.initialize();
  #endif
  ppredc.allocate(1,ntime,1,ntallas,"ppredc");
  #ifndef NO_AD_INITIALIZE
    ppredc.initialize();
  #endif
  pred_CPUE.allocate(1,ntime,"pred_CPUE");
  #ifndef NO_AD_INITIALIZE
    pred_CPUE.initialize();
  #endif
  pred_Bcru.allocate(1,ntime,"pred_Bcru");
  #ifndef NO_AD_INITIALIZE
    pred_Bcru.initialize();
  #endif
  pred_mph.allocate(1,ntime,"pred_mph");
  #ifndef NO_AD_INITIALIZE
    pred_mph.initialize();
  #endif
  Fspr.allocate(1,nedades,"Fspr");
  #ifndef NO_AD_INITIALIZE
    Fspr.initialize();
  #endif
  Zspr.allocate(1,nedades,"Zspr");
  #ifndef NO_AD_INITIALIZE
    Zspr.initialize();
  #endif
  Nspro.allocate(1,nedades,"Nspro");
  #ifndef NO_AD_INITIALIZE
    Nspro.initialize();
  #endif
  Nspr.allocate(1,nedades,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  Nmed.allocate(1,nedades,"Nmed");
  #ifndef NO_AD_INITIALIZE
    Nmed.initialize();
  #endif
  Bspro.allocate("Bspro");
  #ifndef NO_AD_INITIALIZE
  Bspro.initialize();
  #endif
  Bspr.allocate("Bspr");
  #ifndef NO_AD_INITIALIZE
  Bspr.initialize();
  #endif
  ratio_spr.allocate(1,npbr,"ratio_spr");
  #ifndef NO_AD_INITIALIZE
    ratio_spr.initialize();
  #endif
  Fmedian.allocate("Fmedian");
  #ifndef NO_AD_INITIALIZE
  Fmedian.initialize();
  #endif
  Fmed.allocate(1,nedades,"Fmed");
  #ifndef NO_AD_INITIALIZE
    Fmed.initialize();
  #endif
  Zmed.allocate(1,nedades,"Zmed");
  #ifndef NO_AD_INITIALIZE
    Zmed.initialize();
  #endif
  Bsprmed.allocate("Bsprmed");
  #ifndef NO_AD_INITIALIZE
  Bsprmed.initialize();
  #endif
  ratio_Fmed.allocate("ratio_Fmed");
  #ifndef NO_AD_INITIALIZE
  ratio_Fmed.initialize();
  #endif
  Bmed.allocate("Bmed");
  #ifndef NO_AD_INITIALIZE
  Bmed.initialize();
  #endif
  Bo.allocate("Bo");
  #ifndef NO_AD_INITIALIZE
  Bo.initialize();
  #endif
  Brms.allocate(1,npbr,"Brms");
  #ifndef NO_AD_INITIALIZE
    Brms.initialize();
  #endif
  RPRrms.allocate(1,ntime,"RPRrms");
  Frpr.allocate(1,ntime,"Frpr");
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
  penalty.allocate("penalty");
  #ifndef NO_AD_INITIALIZE
  penalty.initialize();
  #endif
  suma4.allocate("suma4");
  #ifndef NO_AD_INITIALIZE
  suma4.initialize();
  #endif
  suma5.allocate("suma5");
  #ifndef NO_AD_INITIALIZE
  suma5.initialize();
  #endif
  suma6.allocate("suma6");
  #ifndef NO_AD_INITIALIZE
  suma6.initialize();
  #endif
  suma7.allocate("suma7");
  #ifndef NO_AD_INITIALIZE
  suma7.initialize();
  #endif
  likeval.allocate(1,9,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  Fpbr.allocate(1,nedades,"Fpbr");
  #ifndef NO_AD_INITIALIZE
    Fpbr.initialize();
  #endif
  Zpbr.allocate(1,nedades,"Zpbr");
  #ifndef NO_AD_INITIALIZE
    Zpbr.initialize();
  #endif
  Np.allocate(1,nedades,"Np");
  #ifndef NO_AD_INITIALIZE
    Np.initialize();
  #endif
  Sp.allocate(1,nedades,"Sp");
  #ifndef NO_AD_INITIALIZE
    Sp.initialize();
  #endif
  Bp.allocate("Bp");
  #ifndef NO_AD_INITIALIZE
  Bp.initialize();
  #endif
  BDp.allocate("BDp");
  #ifndef NO_AD_INITIALIZE
  BDp.initialize();
  #endif
  NMDp.allocate(1,ntallas,"NMDp");
  #ifndef NO_AD_INITIALIZE
    NMDp.initialize();
  #endif
  NV.allocate(1,nedades,"NV");
  #ifndef NO_AD_INITIALIZE
    NV.initialize();
  #endif
  NVtallas.allocate(1,ntallas,"NVtallas");
  #ifndef NO_AD_INITIALIZE
    NVtallas.initialize();
  #endif
  CTP.allocate(1,ntallas,"CTP");
  #ifndef NO_AD_INITIALIZE
    CTP.initialize();
  #endif
  YTP.allocate(1,npbr,"YTP");
  #ifndef NO_AD_INITIALIZE
    YTP.initialize();
  #endif
  CBA.allocate(1,npbr,"CBA");
  NVp.allocate(1,nedades,"NVp");
  #ifndef NO_AD_INITIALIZE
    NVp.initialize();
  #endif
  NVptallas.allocate(1,ntallas,"NVptallas");
  #ifndef NO_AD_INITIALIZE
    NVptallas.initialize();
  #endif
  CTPp.allocate(1,ntallas,"CTPp");
  #ifndef NO_AD_INITIALIZE
    CTPp.initialize();
  #endif
  YTPp.allocate(1,ntime_sim,1,npbr,"YTPp");
  #ifndef NO_AD_INITIALIZE
    YTPp.initialize();
  #endif
  CBAp.allocate(1,npbr,"CBAp");
  SSBp.allocate(1,ntime_sim,1,npbr,"SSBp");
  BTp.allocate(1,ntime_sim,1,npbr,"BTp");
  #ifndef NO_AD_INITIALIZE
    BTp.initialize();
  #endif
  nm1.allocate("nm1");
  #ifndef NO_AD_INITIALIZE
  nm1.initialize();
  #endif
  cuenta1.allocate("cuenta1");
  #ifndef NO_AD_INITIALIZE
  cuenta1.initialize();
  #endif
  nm2.allocate("nm2");
  #ifndef NO_AD_INITIALIZE
  nm2.initialize();
  #endif
  cuenta2.allocate("cuenta2");
  #ifndef NO_AD_INITIALIZE
  cuenta2.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 yrs     = column(data,1);
 Desemb  = column(data,2);
 CPUE    = column(data,3);//cpue
 Bcru    = column(data,4);//reclan
 mph     = column(data,5);//mph
 cv_index(1) = column(data,6);
 cv_index(2) = column(data,7);
 cv_index(3) = column(data,8);
 cv_index(4) = column(data,9);
 nm(1)       = column(data,10);
 nm(2)       = column(data,11);
 dt_C    = column(data,12);
 Frms_bloque=column(data,13);
 Unos_edad   = 1;// lo uso en operaciones matriciales con la edad
 Unos_anos   = 1;// lo uso en operaciones matriciales con el a?o
 Unos_tallas = 1;// lo uso en operaciones matriciales con el a?o
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{5000, 10000, 100000, 500000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-6,1e-7,1e-8, 1e-8}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::userfunction(void)
{
  f =0.0;
 Eval_prob_talla_edad();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_PBR();
 Eval_deinteres();
 Eval_logverosim();
 Eval_funcion_objetivo();
 if(last_phase()){
 	Eval_CTP();
 }
  if(mceval_phase())
    {
    ofstream mcmc_report("mcmc2.mcmc",ios::app);
    mcmc_report << BD(ntime) << " " << BV(ntime) << endl;
    mcmc_report.close();
    ofstream doris("ctp.mcmc",ios::app);
    doris << CBAp(1) << " " << CBAp(2) << " " << CBAp(3) << endl;
    doris.close();
    }
}

void model_parameters::Eval_prob_talla_edad(void)
{
 int i, j;
 Linf = mfexp(log_Linfprior);
 k    = mfexp(log_kprior);
 Lo   = mfexp(log_Lo);
   mu_edad(1) = Lo; 
   for (i=2;i<=nedades;i++){
   mu_edad(i) = mu_edad(i-1)*mfexp(-k)+Linf*(1-mfexp(-k));} //Shnute y Fornier 1980
   sigma_edad = mfexp(log_alfa)+mfexp(log_beta)*mu_edad;
   Prob_talla = ALK( mu_edad, sigma_edad, Tallas);
}

dvar_matrix model_parameters::ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
{
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);
}

void model_parameters::Eval_selectividad(void)
{
 int i,j;
 for (j=1;j<=nbloques1;j++){
   S1(j)=1/(1+exp(-log(19)*(edades-exp(log_L50f(j)))/exp(log_rangof(j))));}
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
            if (yrs(i)>=ybloques1(j)){
                  Sel_f(i)=S1(j);}}}
    Sel_c=1.0;
  if(active(log_L50c)){
  for (j=1;j<=nbloques2;j++){
     S2(j)=1/(1+exp(-log(19)*(edades-exp(log_L50c(j)))/exp(log_rangoc(j))));}
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques2;j++){
           if (yrs(i)>=ybloques2(j)){
                Sel_c(i)=S2(j);}}}
  }
}

void model_parameters::Eval_mortalidades(void)
{
 M=exp(log_M);
 F=elem_prod(Sel_f,outer_prod(mfexp(log_F),Unos_edad));
 Z=F+M;
 S=mfexp(-1.0*Z);
}

void model_parameters::Eval_abundancia(void)
{
 int i, j;
  if(opt_Ro<0)
  {
  log_Ro=log_priorRo;
  }
 No(1) = mfexp(log_Ro); 
 for (int j=2;j<=nedades;j++)
     {
	 No(j)=No(j-1)*mfexp(-1.*M);
	 }
     No(nedades)=No(nedades)/(1-exp(-1.*M));
     SSBo = sum(elem_prod(No*exp(-dt(1)*M)*Prob_talla,elem_prod(msex,Wmed)));
     phi  = SSBo/mfexp(log_Ro);
 h = mfexp(log_h_prior);
 if(opt_Reclu==1) //ByHolt
 {
 alfa = 4*h*mfexp(log_Ro)/(5*h-1);//
 beta = (1-h)*SSBo/(5*h-1);
 }
 if(opt_Reclu==2) // Ricker
 {
 alfa = 1.25*log(5*h)-log(phi);
 beta = 1.25*log(5*h)/SSBo;
 }
 Neq(1) = No(1);
 for (j=2;j<=nedades;j++)
   {
   Neq(j)=Neq(j-1)*mfexp(-Z(1,j-1));
   }
   Neq(nedades)=Neq(nedades)/(1-mfexp(-1.*Z(1,nedades)));
 //---------------------------------------------------------------------------
 // Abundancia inicial
 //---------------------------------------------------------------------------
 N(1)     = elem_prod(Neq,mfexp(dev_log_No));
 BD(1)    = sum(elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1)))*Prob_talla,elem_prod(msex,Wmed)));
 Rpred(1) = N(1,1);
 for (i=1;i<ntime;i++)
 {
     if(opt_Reclu==1) //ByHolt
	 {
     Rpred(i+1)=alfa*BD(i)/(beta+BD(i));
	 }
     if(opt_Reclu==2) //Ricker
	 {
     Rpred(i+1)=BD(i)*mfexp(alfa-1.*beta*BD(i));
	 }
	if(opt_Reclu==3) //sin relaci?n 
     { 
	 Rpred(i+1)= mfexp(log_Ro); 
	 }
     N(i+1,1)          = Rpred(i+1)*mfexp(dev_log_Ro(i+1));
     N(i+1)(2,nedades) = ++elem_prod(N(i)(1,nedades-1),S(i)(1,nedades-1));
     N(i+1,nedades)    = N(i+1,nedades)+N(i,nedades)*S(i,nedades);// grupo plus
     BD(i+1)           = sum(elem_prod(elem_prod(N(i+1),exp(-dt(1)*Z(i+1)))*Prob_talla,elem_prod(msex,Wmed)));
	 //PHI(i)          = BD(i)/N(i+1,1);
 }
     Reclutas          = column(N,1);
     Ntallas           = N*Prob_talla;
}

void model_parameters::Eval_biomasas(void)
{
 NVflo   = elem_prod(elem_prod(N,mfexp(-Z)),Sel_f)*Prob_talla;            //flota
 for(int i=1;i<=ntime;i++){
 NVcru(i)=elem_prod(elem_prod(N(i),mfexp(-dt_C(i)*(Z(i)))),Sel_c(i))*Prob_talla;
 }
 BMflo   = NVflo*Wmed;
 BMcru   = NVcru*Wmed;
 BV      = BMflo;
 BT      = N*Prob_talla*Wmed;  
}

void model_parameters::Eval_capturas_predichas(void)
{
 pred_Ctot_a = elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
 pred_Ctot   = pred_Ctot_a*Prob_talla;
 pred_Desemb = pred_Ctot*Wmed;
 pobs_f      = elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred_f     = elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));
 pobsc       = elem_div(Ncru,outer_prod(rowsum(Ncru+1e-10),Unos_tallas));
 ppredc      = elem_div(NVcru,outer_prod(rowsum(NVcru+1e-10),Unos_tallas));
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
                 pred_Bcru(i)=mfexp(log_qcru(j))*BMcru(i);}
       }
   }
}

void model_parameters::Eval_PBR(void)
{
 if(opt_Ro<0)
 {
 log_Ro=log_priorRo;
 }
 for (int i=1;i<=npbr;i++){
  //Fspr = Sel_f(1)*mfexp(log_Fref(i)); //bloque 1 2002-2009
  //Fspr = Sel_f(9)*mfexp(log_Fref(i)); //bloque 2 2010-2012
    Fspr = Sel_f(ntime)*mfexp(log_Fref(i)); //bloque 3 2013-2020
    Zspr = Fspr+M;
 Nspro(1)=mfexp(log_Ro);
 Nspr(1)=mfexp(log_Ro);
 for (int j=2;j<=nedades;j++)
 {
   Nspro(j)=Nspro(j-1)*mfexp(-1.*M);
   Nspr(j)=Nspr(j-1)*mfexp(-Zspr(j-1));
 }
   Nspro(nedades)=Nspro(nedades)/(1-mfexp(-1.*M)); 
   Nspr(nedades)=Nspr(nedades)/(1-mfexp(-Zspr(nedades))); 
 Bspro   = sum(elem_prod(Nspro*mfexp(-dt(1)*M)*Prob_talla,elem_prod(msex,Wmed)));
 Bspr    = sum(elem_prod(elem_prod(Nspr,mfexp(-dt(1)*Zspr))*Prob_talla,elem_prod(msex,Wmed)));
 ratio_spr(i)=Bspr/Bspro;
 Bo    =  Bspro;
 Brms(i)=Bo*(ratio_spr(i)-0.05);
 }
   RPRrms = BD/(Bo*0.55);
   //Frpr   = elem_div(mfexp(log_F),Frms_bloque);
    Frpr   = exp(log_F)/mfexp((log_Fref(3))); 
}

void model_parameters::Eval_deinteres(void)
{
 Lf_obs   = Tallas*trans(pobs_f);
 Lf_pred  = Tallas*trans(ppred_f);
 Lc_obs   = Tallas*trans(pobsc);
 Lc_pred  = Tallas*trans(ppredc);
 Nv=N;// solo para empezar los calculos
 for (int i=1;i<ntime;i++)
 {
     Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*M);
     Nv(i+1,nedades)=Nv(i+1,nedades)+Nv(i,nedades)*exp(-1.0*M);// grupo plus
 }
 NDv   = elem_prod((Nv*exp(-dt(1)*M))*Prob_talla,outer_prod(Unos_anos,msex));
 BDo   = NDv*Wmed;
 RPR   = elem_div(BD,BDo);
 RPRlp = BD/SSBo;
}

void model_parameters::Eval_logverosim(void)
{
 int i;
 suma1=0; suma2=0;  suma3=0; penalty=0;
 for (i=1;i<=ntime;i++)
 {
   if (CPUE(i)>0){
    suma1+=square(log(CPUE(i)/pred_CPUE(i))*1/cv_index(2,i));}
  if (Bcru(i)>0){
    suma2+=square(log(Bcru(i)/pred_Bcru(i))*1/cv_index(3,i));}
  if (Desemb(i)>0){
    suma3+=square(log(Desemb(i)/pred_Desemb(i))*1/cv_index(1,i));}
  }
}

void model_parameters::Eval_funcion_objetivo(void)
{
 suma4=0; suma5=0; suma6=0; suma7=0; penalty=0; 
 int i;
 likeval(1)=0.5*suma1;//CPUE
 likeval(2)=0.5*suma2;//RECLAN
 likeval(3)=0.5*suma3;// desemb
 suma4 = sum(nm(1)*elem_prod(pobs_f,log(ppred_f)));
 suma5 = sum(nm(2)*elem_prod(pobsc,log(ppredc)));
 likeval(4) = -1.*suma4; //Flota
 likeval(5) = -1.*suma5; //RECLAN
 if(active(dev_log_Ro)){
 likeval(6)=1./(2*square(cvar(1)))*norm2(dev_log_Ro);}
 if(active(dev_log_No)){
 likeval(7)=1./(2*square(cvar(2)))*norm2(dev_log_No);}
 if (active(log_Lo)){
 likeval(8)=0.5*square((log_Loprior-log_Lo)/cv_priors(3));}
 if(active(log_M)){
 penalty+=100*(square(log_M_prior-log_M));}
 if(active(log_qflo)){
 penalty+=0.5*norm2((log_qflo-log_prior_qf)/cv_qf);}
 if(active(log_qcru)){
 penalty+=0.5*norm2((log_qcru-log_prior_qc)/cv_qcru);}
 if(active(log_Fref)){
 penalty+=1000*norm2(ratio_spr-pbr);}
 if (active(log_F)){
 penalty+=1000*(square(log_F(1)-mean(log_F))+square(log_F(2)-mean(log_F))); }
 f=opt_sim*(sum(likeval)+penalty);
}

void model_parameters::Eval_CTP(void)
{
 for (int j=1;j<=npbr;j++){
     Fpbr   = Sel_f(ntime)*mfexp((log_Fref(j)));//
     Zpbr   = Fpbr+M;
	   NV     = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),N(ntime)));
	   NVtallas=NV*Prob_talla;
     CTP    = elem_prod(NVtallas,Wmed);
     YTP(j) = sum(CTP);
     }
 for (int j=1;j<=npbr;j++){
 Np   = N(ntime);
 Sp   = S(ntime);
 BDp  = BD(ntime);
 Fpbr = F(ntime);//
 Zpbr = Z(ntime);
  if(opt_Ro<0)
  {
  log_Ro=log_priorRo;
  }
   for (int i=1;i<=ntime_sim;i++)
   {
  Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
  Np(nedades)+=Np(nedades)*Sp(nedades);
   if(oprec==1){
      Np(1)=(mfexp(log_Ro));} // Estima CTP con R_alto (2012)
   if(oprec==2){
      Np(1)=(Reclutas(11));} // Estima CTP con R_alto (2012)
   if(oprec==3){
      Np(1)=(Reclutas(17));}     // Estima CTP con R_bajo (2018)
   if(oprec==4){                // reclutamiento para proyecci?n seg?n relaci?n stock-recluta o no
      if(opt_Reclu==1){ //Beverton Y Holt
       Np(1)=pR*alfa*BDp/(beta+BDp)*mfexp(0.5*square(cvar(1)));}
      if(opt_Reclu==2){ //Ricker
       Np(1)=pR*BDp*mfexp(alfa-1.0*beta*BDp)*mfexp(0.5*square(cvar(1)));}
      if(opt_Reclu==3){ //sin relaci?n stock-recluta
       Np(1)=pR*mfexp(log_Ro+0.5*square(cvar(1)));}}	  
  Fpbr = Sel_f(ntime)*mfexp((log_Fref(j)));//
  Zpbr = Fpbr+M;
  Sp=exp(-1.*Zpbr);
  NVp  = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(Np,(1-Sp)));
  NVptallas=NVp*Prob_talla;
  CTPp = elem_prod(NVptallas,Wmed);
  YTPp(i,j)=sum(CTPp);
  Bp=sum(elem_prod(Np*Prob_talla,Wmed));
  NMDp = elem_prod(Np,mfexp(-dt(1)*(Zpbr)))*Prob_talla;
  BDp  = sum(elem_prod(elem_prod(NMDp,msex),Wmed));
  SSBp(i,j)=BDp;
  BTp(i,j)=Bp;
 }
  CBA(j) = YTP(j);
  CBAp(j) = YTPp(1,j);// es para el a?o proyectado revisar si es 1 o 2 
 }
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
 report << "YRS" << endl;
 report << yrs << endl;
 report << "desemb" << endl;
 report << Desemb << endl;
 report << "desemb_pred" << endl;
 report << pred_Desemb << endl;
 report << "reclan" << endl;
 report << Bcru << endl;
 report << "reclan_pred" << endl;
 report << pred_Bcru << endl;
 report << "cpue" << endl;
 report << CPUE << endl;
 report << "cpue_pred" << endl;
 report << pred_CPUE << endl;
 report << "mph" << endl;
 report << mph << endl;
 report << "mph_pred" << endl;
 report << pred_mph << endl;
 report << "BD" << endl;
 report << BD << endl;
 report << "BT" << endl;
 report << BT << endl;
 report << "BV" << endl;
 report << BMflo << endl;
 report << "Rech_pre_est" << endl;
 report << Rpred<< endl;
 report << "Reclutas" << endl;
 report << column(N,1)<< endl;
 report << "F" << endl;
 report << mfexp(log_F)<< endl;
 report << "Lf_obs" << endl;
 report << Lf_obs<< endl;
 report << "Lf_pred" << endl;
 report << Lf_pred<< endl;
 report << "Lc_obs" << endl;
 report << Lc_obs<< endl;
 report << "Lc_pred" << endl;
 report << Lc_pred<< endl;
 //-------------------------------------------
 report << "Sflo_age" <<endl;
 report << Sel_f << endl;
 report << "Scru_age" << endl;
 report << Sel_c << endl;
 //-------------------------------------------
 report << "pf_obs" << endl;
 report << (pobs_f)<< endl;
 report << "pf_pred" << endl;
 report << (ppred_f)<< endl;
 report << "pobs_RECLAN" << endl;
 report << (pobsc)<< endl;
 report << "ppred_RECLAN" << endl;
 report << (ppredc)<< endl;
 //------------------------------------------
 report << "N" << endl;
 report << N << endl;
 report << "Ntallas" << endl;
 report << Ntallas << endl;
 report << "Capt_age" <<endl;
 report << pred_Ctot <<endl;
 report << "Capt_talla" <<endl;
 report << pred_Desemb <<endl;
 report << "F_age" <<endl;
 report << F << endl;
 report << "BDo" << endl;
 report << BDo << endl;
 report << "BDoLP" << endl;
 report << SSBo << endl;
 report << "RPR" << endl;
 report << RPR << endl;
 report << "BD/BDo" << endl;
 report << RPRlp << endl;
 report << "RPRrms"<<endl;
 report <<  RPRrms << endl;
 report << "Frpr"<<endl;
 report <<  Frpr <<endl;
 report << "log_Ro" <<endl;
 report << log_Ro <<endl;
 report << "Lo"  << endl;
 report << exp(log_Lo) << endl;
 report << "M" << endl;
 report << M <<endl;
 report << "mu_edad" << endl;
 report << mu_edad << endl;
 report << "Prob_talla" << endl;
 report << Prob_talla << endl;
 report << "likeval cpue Crucero   Desemb   propflo  pobsc  Ro    No  Lo" <<endl;
 report << likeval << endl;
 report << "q_flo" <<endl;
 report << exp(log_qflo) << endl;
 report << "q_cru_reclan" <<endl;
 report << exp(log_qcru) << endl;
 report << "alfa" << endl;
 report << alfa << endl;
 report << "beta" <<endl;
 report << beta << endl;
 report << "Bo"<<endl;
 report << Bo << endl;
 report << "Brms" << endl;
 report << Brms << endl;
 report << "Rproy" << endl;
 report << Np(1) << endl;
 report << "BT_proy" << endl;
 report << BTp << endl;
 report << "BD_proy" << endl;
 report << SSBp << endl;
 //CBA actual
 report << "NV" << endl;
 report << NV << endl;
 report << "NVtallas" << endl;
 report << NVtallas << endl;
 report << "CTP" << endl;
 report << CTP << endl;
 report << "YTP" << endl;
 report << YTP << endl;
 report << "CBA_actual" << endl;
 report << CBA << endl;
  //CBA proyectada
 report << "N_proyect" << endl; // abundancia del ultimo a?o
 report << Np << endl;
 report << "NVp" << endl;
 report << NVp << endl;
 report << "NVptallas" << endl;
 report << NVptallas << endl;
 report << "CTPp" << endl;
 report << CTPp << endl;
 report << "YTPp" << endl;
 report << YTPp << endl;
 report << "C_proy" << endl;
 report << CBAp << endl;
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
 //######################################################################
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
 arrmblsize = 50000000;
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
