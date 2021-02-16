%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates structured coalescent xmls from the master trees. Always creates
% 3 xmls per tree with different initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% create lisco files
tree_files = dir('master/*.tree');

system('rm -r xmls');
system('mkdir xmls');

states = 3;

for i = 1 : length(tree_files)
    disp(tree_files(i).name)
    % read tree files
    g = fopen(['master/' tree_files(i).name],'r');
    t = textscan(g,'%s'); fclose(g);
    
   % coalescing
    tree_tmp2 = regexprep(t{1}(end-1),'&type="L",location="(\d*)",reaction="Coalescence",time=(\d*).(\d*)','');

    %migrating
    tree_tmp1 = regexprep(tree_tmp2,'&type="L",location="(\d*)",reaction="Migration",time=(\d*).(\d*)','');

    
    % make the MASTER tree compatible with BEAST2
    % sampling
    tree_tmp1 = regexprep(tree_tmp1,'E[-](\d)]',']');
    tip_locs = regexp(tree_tmp1,'[&type="L",location="(\d*)",time=(\d*)\.(\d*)\]','match');
     
    for j = 1 : length(tip_locs{1})
        loc = regexprep(tip_locs{1}{j},'[&type="L",location="(\d*)",time=(\d*)\.(\d*)\]','$1');
        tree_tmp1 = strrep(tree_tmp1,tip_locs{1}{j},['loc_' loc 'kickout']);
        tree_tmp1 = strrep(tree_tmp1,'kickout','');
    end

    tree_tmp = regexprep(tree_tmp1,'(\d*)loc_','inv$1loc_');
    
    tree = strrep(tree_tmp{1},'[]','');
    if ~isempty(strfind(tree,']'))
        b = strfind(tree,']');
        c = tree((b-50):(b+50));
        disp(tree_files(i).name)
    end

    % get the leafnames
    ptree = phytreeread(tree);
    leafnames = get(ptree,'leafnames');
    
    print_tree = tree;
    
    % get the covariates
    
    % make tripletts of all runs with different random initial values
    for tr = 1 : 1    
        % make the xmls for the structcoal
        flog = strrep(tree_files(i).name,'master.tree',sprintf('%dmascot',tr));
        fname = sprintf('xmls/%s.xml',flog);
        f = fopen(fname,'w');


        fprintf(g,'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate=''Standard'' beautistatus='''' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.mascot.dynamics:beast.mascot.distribution:beast.mascot.logger" version="2.0">\n');

        fprintf(g,'\t<data id="sequences" name="alignment">\n');
        for j = 1 : length(leafnames)
            fprintf(g,'\t\t<sequence id="%s" taxon="%s" totalcount="4" value="??"/>\n',leafnames{j},leafnames{j});
        end
        fprintf(g,'\t</data>\n');
        fprintf(g,'<map name="prior">beast.math.distributions.Prior</map>\n');
        fprintf(g,'<map name="Normal">beast.math.distributions.Normal</map>\n');

        fprintf(g,'\t<run id="mcmc" spec="MCMC" chainLength="250000">\n');
        fprintf(g,'\t\t<state id="state" storeEvery="5000">\n');
        fprintf(g,'\t\t\t<stateNode id="tree" spec="beast.app.mascot.beauti.TreeWithTrait">\n');
        fprintf(g,'\t\t\t\t<typeTrait id="typeTraitSet.t" spec="beast.evolution.tree.TraitSet" traitname="type" value="');
        for j = 1 : length(leafnames)-1
            tmp1 = strsplit(leafnames{j},'_');
            fprintf(g,'%s=state%s,',leafnames{j},tmp1{end});
        end
        tmp1 = strsplit(leafnames{end},'_');
        fprintf(g,'%s=state%s">\n',leafnames{end},tmp1{end});

        fprintf(g,'\t\t\t\t\t<taxa id="TaxonSet.0" spec="TaxonSet">\n');
        fprintf(g,'\t\t\t\t\t\t<alignment idref="sequences"/>\n');
        fprintf(g,'\t\t\t\t\t</taxa>\n');
        fprintf(g,'\t\t\t\t</typeTrait>\n');
        fprintf(g,'\t\t\t</stateNode>\n');        
        
        fprintf(g,'\t\t\t<parameter id="Nenull.1" name="stateNode" dimension="1">0</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="Nenull.2" name="stateNode" dimension="1">0</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="Nenull.3" name="stateNode" dimension="1">0</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="growth.1" name="stateNode" dimension="1">0.2</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="growth.2" name="stateNode" dimension="1">0.1</parameter>\n');
        fprintf(g,'\t\t\t<parameter id="growth.3" name="stateNode" dimension="1">0</parameter>\n');

        
        fprintf(g,'\t\t\t<parameter id="migration" name="stateNode" dimension="0">0.1</parameter>\n');
        fprintf(g,'\t\t</state>\n');
        

                
        fprintf(g,'\t\t<init spec="beast.util.TreeParser" id="NewickTree.t:Species" adjustTipHeights="false"\n');
        fprintf(g,'\t\t\tinitial="@tree" taxa="@sequences"\n');
        fprintf(g,'\t\t\tIsLabelledNewick="true"\n');
        fprintf(g,'\t\t\tnewick="%s"/>\n',print_tree);
        fprintf(g,' \t\t<distribution id="posterior" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t<distribution id="prior" spec="util.CompoundDistribution">\n');
        
%         fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@NeMean">\n');
%         fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.OneOnX"/>\n');
%         fprintf(g,'\t\t\t\t</distribution>\n');

        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@migration">\n');
        fprintf(g,'\t\t\t\t<distr spec="beast.math.distributions.Normal"   mean="-1" sigma="1.0"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');


        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@Nenull.1">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0" sigma="1.0"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');

        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@Nenull.2">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0" sigma="1.0"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');

        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@Nenull.3">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0" sigma="1.0"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');



        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@growth.1">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0.5" sigma="0.25"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');

        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@growth.2">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0.5" sigma="0.25"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');

        fprintf(g,'\t\t\t\t<distribution spec=''beast.math.distributions.Prior'' x="@growth.3">\n');
        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0.5" sigma="0.25"/>\n');
        fprintf(g,'\t\t\t\t</distribution>\n');


                
       

        fprintf(g,'\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t<distribution id="likelihood" spec="util.CompoundDistribution">\n');
        fprintf(g,'\t\t\t\t<distribution id="coalescent" spec="Mascot">\n');
        fprintf(g,'\t\t\t\t\t<structuredTreeIntervals spec="StructuredTreeIntervals" id="TreeIntervals" tree="@tree"/>\n');        
        fprintf(g,'\t\t\t\t\t<dynamics id="piecewiseConstant" spec="beast.mascotskyline.dynamics.DynamicEffectivePopulationSizes" forwardsMigration="@migration" dimension="0" typeTrait="@typeTraitSet.t">\n');
        fprintf(g,'\t\t\t\t\t\t<NeDynamics id="expo.1" spec="beast.mascotskyline.parametric.Exponential" NeNull="@Nenull.1" growthRate="@growth.1"/>\n');
        fprintf(g,'\t\t\t\t\t\t<NeDynamics id="expo.2" spec="beast.mascotskyline.parametric.Exponential" NeNull="@Nenull.2" growthRate="@growth.2"/>\n');
        fprintf(g,'\t\t\t\t\t\t<NeDynamics id="expo.3" spec="beast.mascotskyline.parametric.Exponential" NeNull="@Nenull.3" growthRate="@growth.3"/>\n');
        fprintf(g,'\t\t\t\t\t\t<rateShifts id="rateshifts.t" spec="beast.mascotskyline.skyline.RateShifts" tree="@tree">\n');
        fprintf(g,'\t\t\t\t\t\t\t<parameter id="rateshifts" name="relativeRateShifts">0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00 1.01</parameter>\n');
        fprintf(g,'\t\t\t\t\t\t</rateShifts>\n');
        fprintf(g,'\t\t\t\t\t</dynamics>\n');        
        fprintf(g,'\t\t\t\t</distribution>\n');
        fprintf(g,'\t\t\t</distribution>\n');
        fprintf(g,'\t\t</distribution>\n');

        
        fprintf(g,'\t\t<operator id="NeCarryinglScaler.1" spec="beast.mascotskyline.operators.AdaptiveMultivariateGaussianOperator" windowSize="0.1" weight="3.0">\n');
        fprintf(g,'\t\t\t<parameter idref="Nenull.1"/>\n');
        fprintf(g,'\t\t\t<parameter idref="growth.1"/>\n');
        fprintf(g,'\t\t\t<parameter idref="Nenull.2"/>\n');
        fprintf(g,'\t\t\t<parameter idref="growth.2"/>\n');
        fprintf(g,'\t\t\t<parameter idref="Nenull.3"/>\n');
        fprintf(g,'\t\t\t<parameter idref="growth.3"/>\n');
        fprintf(g,'\t\t\t<parameter idref="migration"/>\n');
%         fprintf(g,'\t\t\t<logTransform id="logTransform" estimate="false" spec="parameter.BooleanParameter">false false false false false false false false false true</logTransform>\n');
        fprintf(g,'\t\t</operator>\n');

        fprintf(g,'\t\t<logger id="probtreelog" fileName="%s.trees" logEvery="5000" mode="tree">\n',flog);
        fprintf(g,'\t\t\t<log id="logTrees" spec="StructuredTreeLogger" mascot="@coalescent"/>\n');
        fprintf(g,'\t\t</logger>\n');

        fprintf(g,'\t\t<logger id="tracelog" fileName="%s.log" logEvery="200" model="@posterior" sanitiseHeaders="true" sort="smart">\n',flog);
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t\t<log idref="prior"/>\n');

        fprintf(g,'\t\t\t<log spec="RootStateLogger" id="RootStateLogger" mascot="@coalescent"/>\n');
        fprintf(g,'\t\t\t<log idref="piecewiseConstant"/>\n');
        fprintf(g,'\t\t\t<log idref="Nenull.1"/>\n');
        fprintf(g,'\t\t\t<log idref="Nenull.2"/>\n');
        fprintf(g,'\t\t\t<log idref="Nenull.3"/>\n');
        fprintf(g,'\t\t\t<log idref="growth.1"/>\n');
        fprintf(g,'\t\t\t<log idref="growth.2"/>\n');
        fprintf(g,'\t\t\t<log idref="growth.3"/>\n');
        fprintf(g,'\t\t</logger>\n');
        
        fprintf(g,'\t\t<logger id="screenlog" logEvery="1000">\n');
        fprintf(g,'\t\t\t<log idref="posterior"/>\n');
        fprintf(g,'\t\t</logger>\n');
        fprintf(g,'\t</run>\n');
        fprintf(g,'</beast>\n');
        fclose(f);
    end
end
