#include "stdafx.h"
#include "kronecker.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("Kronecker graphs. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  Env = TEnv(argc, argv, TNotify::StdNotify);
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../as20graph.txt", "Input graph file (single directed edge per line)");
  const TInt NZero = Env.GetIfArgPrefixInt("-n0:", 2, "Matrix size");
  const TStr InitMtx = Env.GetIfArgPrefixStr("-m:", "0.9 0.7; 0.5 0.2", "Kronecker matrix").GetLc();
  const TInt WarmUp =  Env.GetIfArgPrefixInt("-w:", 10000, "Samples to warm up");
  const TInt NSamples = Env.GetIfArgPrefixInt("-s:", 100000, "Samples per log-liklihood estimation");
  const TStr Perm = Env.GetIfArgPrefixStr("-p:", "d", "Initial node permutation: d:Degree, r:Random, o:Order").GetLc();
  //const TInt GradType = Env.GetIfArgPrefixInt("-gt:", 1, "1:Grad1, 2:Grad2");
  const bool ScaleInitMtx = Env.GetIfArgPrefixBool("-sim:", false, "Scale the initiator to match the number of edges");
  const TFlt PermSwapNodeProb = Env.GetIfArgPrefixFlt("-nsp:", 1.0, "Probability of using NodeSwap (vs. EdgeSwap) MCMC proposal distribution");
  
  // load graph
  PNGraph G;
  if (InFNm.GetFExt().GetLc()==".ungraph") {
    TFIn FIn(InFNm);  G=TSnap::ConvertGraph<PNGraph>(TUNGraph::Load(FIn), true); }
  else if (InFNm.GetFExt().GetLc()==".ngraph") {
    TFIn FIn(InFNm);  G=TNGraph::Load(FIn); }
  else {
    G = TSnap::LoadEdgeList<PNGraph>(InFNm, 0, 1);
  }
  
  // calculate liklihood
  TKronMtx InitKronMtx = TKronMtx::GetMtx(InitMtx);
  InitKronMtx.Dump("KRONECKER PARMETERS", false);
  TKroneckerLL KronLL(G, InitKronMtx, PermSwapNodeProb);
  
  if (ScaleInitMtx) {
    InitKronMtx.SetForEdges(G->GetNodes(), G->GetEdges()); }
  KronLL.InitLL(G, InitKronMtx);
  InitKronMtx.Dump("SCALED KRONECKER PARMETERS", false);
  
  KronLL.SetPerm(Perm.GetCh(0));
  double LogLike = 0;
  TFltV CurGradV;
  
  KronLL.SampleGradient(WarmUp, NSamples, LogLike, CurGradV);

  fprintf(stdout, "Input\t%s\n", InFNm.CStr());
  TStrV ParamV; Env.GetCmLn().SplitOnAllCh(' ', ParamV);
  fprintf(stdout, "Command line options\n");
  for (int i = 0; i < ParamV.Len(); i++) {
    fprintf(stdout, "\t%s\n", ParamV[i].CStr()+(ParamV[i][0]=='-'?1:0)); }
  fprintf(stdout, "Loglikelihood\t%18.16e\n", LogLike);
  fprintf(stdout, "Gradient:\n");
  for (int p = 0; p < KronLL.GetParams(); p++) {
      fprintf(stdout, "\t%18.16e\n", (double)CurGradV[p]); }
  fprintf(stdout, "RunTime\t%g\n", ExeTm.GetSecs());

  Catch
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  return 0;
}
