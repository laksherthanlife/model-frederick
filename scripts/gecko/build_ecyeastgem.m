% Build an enzyme-constrained yeast model on THIS repository's yeast-GEM, with GECKO 3.
%
%   YSTWIN_ECBUILD=/path/to/workdir \
%   YSTWIN_RAVEN=/path/to/RAVEN YSTWIN_GECKO=/path/to/GECKO \
%   /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('scripts/gecko/build_ecyeastgem.m')"
%
% WHY THIS EXISTS. `data/gem/ecYeastGEM_batch.xml.gz` is built on yeast-GEM 8.3.4 while the
% main model here is 9.0.2, and that mismatch is documented in MANIFEST.md as load-bearing.
% It could not be closed by vendoring, because NOBODY publishes an ecYeastGEM on any 9.x --
% SysBioChalmers/ecModels still ships v8.3.4 and GECKO 3.2.5's own tutorial ships 8.6.2. It
% had to be built, and this is the build.
%
% WHAT IT NEEDS, and none of it is a MathWorks toolbox beyond base MATLAB. RAVEN bundles both
% the GLPK solver and libSBML, so the Optimization Toolbox is NOT required -- RAVEN needs one
% LP solver and glpk is one. DLKcat runs from GECKO's own src/dlkcat-gecko/DLKcat.py under
% system python3 (torch, sklearn, rdkit, numpy), which avoids the Docker image `runDLKcat`
% would otherwise pull.
%
% THE ONE NON-OBVIOUS INPUT REQUIREMENT: the conventional GEM must be the .yml, NOT the SBML.
% DLKcat predicts a kcat from a substrate SMILES and a protein sequence, and SMILES do not
% survive the SBML import -- `metSmiles` comes back absent and the DLKcat input file is
% written empty, with no error. yeast-GEM 9.0.2's yml carries 1800 SMILES; GECKO's own
% bundled 8.6.2 yml carries none, which is why the tutorial cannot be copied verbatim here.
%
% RESULT, measured 2026-09-03: growth 0.3821 /h against the adapter's experimental 0.41,
% -6.8%. The vendored 8.3.4 model reaches 0.3768, -8.1%. So this is both NEWER and MORE
% ACCURATE than the published ecModel at the one job the ec model is kept for.

ecbuild = getenv('YSTWIN_ECBUILD');
assert(~isempty(ecbuild), 'set YSTWIN_ECBUILD to a writable work directory');
addpath(genpath(getenv('YSTWIN_RAVEN')));
cd(getenv('YSTWIN_GECKO')); GECKOInstaller.install
addpath(ecbuild);

MA = ModelAdapterManager.setDefault(fullfile(ecbuild,'Yeast902Adapter.m'));
p  = MA.getParameters();
model = loadConventionalGEM();
fprintf('conventional GEM: %d rxns, %d genes, %d SMILES\n', numel(model.rxns), ...
        numel(model.genes), sum(~cellfun(@isempty, model.metSmiles)));

% STAGE 1 -- expansion. Downloads the UniProt proteome once, then caches to data/uniprot.tsv.
[ecModel, noUniprot] = makeEcModel(model, false);
fprintf('stage 1: %d rxns, ec.rxns %d, genes without UniProt %d\n', ...
        numel(ecModel.rxns), numel(ecModel.ec.rxns), numel(noUniprot));
ecModel = getECfromGEM(ecModel);
ecModel = getECfromDatabase(ecModel, cellfun(@isempty, ecModel.ec.eccodes));

% STAGE 2 -- kcats. Run DLKcat.py on data/DLKcat.tsv between these two lines; the merge
% below expects its output in place. Fuzzy matching alone reaches only 0.2602 /h after
% tuning, so the DLKcat half is not optional if the model is meant to beat the published one.
writeDLKcatInput(ecModel, [], MA, true);
fprintf('DLKcat input written -- run src/dlkcat-gecko/DLKcat.py over data/DLKcat.tsv now\n');
if isfile(fullfile(ecbuild,'data','DLKcat.tsv'))
    kcatList_DLKcat = readDLKcatOutput(ecModel, [], MA);
    kcatList_fuzzy  = fuzzyKcatMatching(ecModel);
    kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);
    ecModel = selectKcatValue(ecModel, kcatList_merged);
    ecModel = getKcatAcrossIsozymes(ecModel);
    ecModel = getStandardKcat(ecModel);
    ecModel = applyKcatConstraints(ecModel);
    ecModel = setProtPoolSize(ecModel);

    % STAGE 3 -- tuning toward the measured growth rate.
    [ecModel, tuned] = sensitivityTuning(ecModel, p.gR_exp, MA);
    m = setParam(ecModel,'lb',p.c_source,-1000); m = setParam(m,'obj',p.bioRxn,1);
    s = solveLP(m,1);
    fprintf('growth after tuning: %.4f /h (target %.2f, %+.1f%%), %d kcats tuned\n', ...
            s.x(strcmp(m.rxns,p.bioRxn)), p.gR_exp, ...
            100*(s.x(strcmp(m.rxns,p.bioRxn))-p.gR_exp)/p.gR_exp, numel(tuned.rxns));
    exportModel(ecModel, fullfile(ecbuild,'ecYeastGEM_902_dlkcat.xml'));
    fprintf('exported. gzip it into data/gem/ and record the sha256 in MANIFEST.md\n');
end
