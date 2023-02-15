%% Pipeline Builder
%
%   Level 1: load hypergraph data
%       user functions: load_HG
%   Level 2: hypergraph observation methods
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: February 9, 2023

clear all; close all;

import bioinfo.pipeline.Pipeline
import bioinfo.pipeline.blocks.*

% Create Pipeline
HGOuniformER = Pipeline;

% Level 1: Load hypergraph
loadHG = UserFunction(@ load_HG,RequiredArguments="ds",OutputArguments="HG");                               % Block 1
addBlock(HGOuniformER,loadHG);

% Level 2: Observe system
OVH = UserFunction(@ observeV, RequiredArguments=["HG","rH","mp"],OutputArguments="OVH");
OVN = UserFunction(@ observeV, RequiredArguments=["HG","rN","mp"],OutputArguments="OVN");
OVL = UserFunction(@ observeV, RequiredArguments=["HG","rL","mp"],OutputArguments="OVL");
OEH = UserFunction(@ observeE, RequiredArguments=["HG","rH","mp"],OutputArguments="OEH");
OEN = UserFunction(@ observeE, RequiredArguments=["HG","rN","mp"],OutputArguments="OEN");
OEL = UserFunction(@ observeE, RequiredArguments=["HG","rL","mp"],OutputArguments="OEL");
addBlock(HGOuniformER,[OVH, OVN, OVL, OEH, OEN, OEL]);

% Connections: Level 1 -- Level 2
connect(HGOuniformER,loadHG,OVH,["HG","HG"]);
connect(HGOuniformER,loadHG,OVN,["HG","HG"]);
connect(HGOuniformER,loadHG,OVL,["HG","HG"]);
connect(HGOuniformER,loadHG,OEH,["HG","HG"]);
connect(HGOuniformER,loadHG,OEN,["HG","HG"]);
connect(HGOuniformER,loadHG,OEL,["HG","HG"]);

% Level 3
collect = UserFunction(@collector, RequiredArguments=["ds", "OVH","OVN","OVL","OEH","OEN","OEL"],OutputArguments="O");
addBlock(HGOuniformER, collect);

% Connections: Level 2 -- Level 3
connect(HGOuniformER,OVH,collect,["OVH","OVH"]);
connect(HGOuniformER,OVN,collect,["OVN","OVN"]);
connect(HGOuniformER,OVL,collect,["OVL","OVL"]);
connect(HGOuniformER,OEH,collect,["OEH","OEH"]);
connect(HGOuniformER,OEN,collect,["OEN","OEN"]);
connect(HGOuniformER,OEL,collect,["OEL","OEL"]);



biopipelineDesigner(HGOuniformER)