% Test inference code

Opts = pilot_test_inference_opts_general();

fn = 'img_0156'; % Load a test image



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load features from different modules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- i.e. Stuff unary features, detector, planes

% 1.1. Stuff Unary features
Feat.StuffUnary = [];
SegFn = fullfile(Opts.FeatStuffUnaryPath, [fn '.seg']);
Seg  = ReadMatrixFromTxtFile(SegFn) + 1;
nSeg = max(Seg(:));

LabelFn = fullfile(Opts.FeatStuffUnaryPath, [fn '.labels.txt']);
Label = ReadMatrixFromTxtFile(LabelFn) + 1;
Label(Label>=14) = 14; % Don't distinguish walls with different orientations

ConfFn = fullfile(Opts.FeatStuffUnaryPath, [fn '.boosted.txt']);
Conf = ReadMatrixFromTxtFile(ConfFn);
Feat.StuffUnary = Conf(:, 2:end)';
Feat.StuffUnary = Feat.StuffUnary - min(Feat.StuffUnary(:)) + 0.01;
Feat.StuffUnary = Feat.StuffUnary / max(Feat.StuffUnary(:));
Feat.StuffUnary = -log(Feat.StuffUnary);
Feat.StuffUnary = Feat.StuffUnary * Opts.w_stuff_unary;

% 1.2. Detection features (X,Y)
% hop_det(i).ind
%           .w
%           .gamma = #superpixel
%           .Q = R * sum(w_i)
%           .thing_uw = score * #superpixel
Feat.HopDet = struct('ind', {}, 'w', {}, 'gamma', {}, 'Q', {}, 'thing_uw', {});
ValidDetIndex = [0 0]; % meaningless set
for iC = 1:numel(Opts.ObjClsList)
    if isfield(Opts, 'DetScoreBiasPerClass'),
        DetScoreBias_backup = Opts.DetScoreBias;
        Opts.DetScoreBias = Opts.DetScoreBiasPerClass{iC};
    end
    if isfield(Opts, 'w_gamma_per_class'),
        w_gamma_backup = Opts.w_gamma;
        Opts.w_gamma = Opts.w_gamma_per_class{iC};
    end
    if isfield(Opts, 'w_thing_uw_per_class'),
        w_thing_uw_backup = Opts.w_thing_uw;
        Opts.w_thing_uw = Opts.w_thing_uw_per_class{iC};
    end
    tDet = load(fullfile(Opts.FeatDetPath, Opts.ObjClsList{iC}, [fn '.mat']));
    for iB = 1:numel(tDet.Bboxes)
        [~, ClassId] = ismember(Opts.ObjClsList{iC}, Opts.AllClsList);
        Feat.HopDet(end+1) = GenFeatHopDet(...
            tDet.Bboxes(iB).det, Seg, ClassId, Opts);
        ValidDetIndex = [ValidDetIndex; iC iB];
        if Feat.HopDet(end).Q == 0 || Feat.HopDet(end).thing_uw < 0,
            Feat.HopDet = Feat.HopDet(1:end-1);
            ValidDetIndex = ValidDetIndex(1:end-1, :);
        end
    end
    if isfield(Opts, 'DetScoreBiasPerClass'),
        Opts.DetScoreBias = DetScoreBias_backup;
    end
    if isfield(Opts, 'w_gamma_per_class'),
        Opts.w_gamma = w_gamma_backup;
    end
    if isfield(Opts, 'w_thing_uw_per_class'),
        Opts.w_thing_uw = w_thing_uw_backup;
    end
    
end
ValidDetIndex = ValidDetIndex(2:end, :);

% 1.3. Plane features (X,P)
Feat.HopPlane = struct('ind', {}, 'w', {}, 'gamma', {}, 'Q', {});
if ~exist(fullfile(Opts.FeatPlanePath, [fn '.mat']), 'file'),
    continue;
end
load(fullfile(Opts.FeatPlanePath, [fn '.mat']), 'Planes');
ValidPlaneIndex = [];
for iP = 1:numel(Planes),
    tP = Planes(iP);
    
    Feat.HopPlane(end+1) = GenFeatHopPlane(tP, Seg, Opts);
    if Feat.HopPlane(end).Q == 0,
        Feat.HopPlane = Feat.HopPlane(1:end-1);
    else
        ValidPlaneIndex = [ValidPlaneIndex iP];
    end
end

% 1.4. Plane x Supp (Thing) features (P,S)
Feat.ThingPlane = struct('thing_id', {}, 'plane_id', {}, 'D', {}, 'n', {}, 'dr', {});
iPS = 0;
for iP = ValidPlaneIndex, %1:numel(Planes),
    tP = Planes(iP);
    iThing = 0;
    for iC = 1:numel(Opts.ObjClsList),
        tDet = load(fullfile(Opts.FeatDetPath, Opts.ObjClsList{iC}, [fn '.mat']));
        for iB = 1:numel(tDet.Bboxes),
			if isempty(find(ValidDetIndex(:,1)==iC & ValidDetIndex(:,2)==iB)),
				continue;
			end
            iPS = iPS + 1;
            iThing = iThing + 1;
            Feat.ThingPlane(iPS).thing_id = iThing;
            Feat.ThingPlane(iPS).plane_id = iP;
            
            
            lo = tDet.Bboxes(iB).cuboid;
            lp = tP.points_3d;
            x1 = lp(1,:); x2 = lp(2,:); x3 = lp(3,:);
            t_n = get_plane_norm(x1,x2,x3);
            t_n = t_n / sqrt(sum(t_n.^2));
            t_op = [(lo(1)+lo(4))/2 (lo(2)+lo(5))/2 lo(3)] - mean(lp);
            t_op = t_op / sqrt(sum(t_op.^2));
            
            lo = tDet.Bboxes(iB).cuboid;
            Feat.ThingPlane(iPS).D = dist([(lo(1)+lo(4))/2 (lo(2)+lo(5))/2 lo(3)], ... %location of object
                                          [mean(lp)]');
            Feat.ThingPlane(iPS).n = t_n .* t_op;
            Feat.ThingPlane(iPS).dr = 1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Use them to construct a graph, and then infer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Feat.Sc = ones(size(Feat.StuffUnary, 1)) - eye(size(Feat.StuffUnary, 1));
Feat.Ss = sparse(zeros(size(Feat.StuffUnary, 2)));

EmptyHopDet = struct('ind', {}, 'w', {}, 'gamma', {}, 'Q', {}, 'thing_uw', {});
EmptyHopPlane = struct('ind', {}, 'w', {}, 'gamma', {}, 'Q', {});
EmptyThingPlane = struct('thing_id', {}, 'plane_id', {}, 'D', {}, 'n', {}, 'dr', {});


% Crf model
disp('Starting Crf model');
gch_crf = GraphCut('open_hop', Feat.StuffUnary, Feat.Sc, Feat.Ss, struct([]), ...
                  EmptyHopDet, struct([]), struct([]), EmptyHopPlane, EmptyThingPlane); 
[gch_crf var_assig_crf.class var_assig_crf.det_label var_assig_crf.plane_label] = GraphCut('expand',gch_crf);
var_assig_crf.class = var_assig_crf.class + 1;
gch_crf = GraphCut('close', gch_crf);
acc_crf =  sum(var_assig_crf.class == Label) / numel(Label);

% Crf+D model
disp('Starting Crf+D model');
gch_crf_d = GraphCut('open_hop', Feat.StuffUnary, Feat.Sc, Feat.Ss, struct([]), ...
                  Feat.HopDet, struct([]), struct([]), EmptyHopPlane, EmptyThingPlane); 
[gch_crf_d var_assig_crf_d.class var_assig_crf_d.det_label var_assig_crf_d.plane_label] = GraphCut('expand',gch_crf_d);
var_assig_crf_d.class = var_assig_crf_d.class + 1;
gch_crf_d = GraphCut('close', gch_crf_d);
acc_crf_d = sum(var_assig_crf_d.class == Label) / numel(Label);

% Crf+D+PX model
disp('==========================================');
disp('Starting Crf+D+PX model');
disp('==========================================');
gch_crf_d_px = GraphCut('open_hop', Feat.StuffUnary, Feat.Sc, Feat.Ss, struct([]), ...
	Feat.HopDet, struct([]), struct([]), Feat.HopPlane, EmptyThingPlane);
[gch_crf_d_px var_assig_crf_d_px.class var_assig_crf_d_px.det_label var_assig_crf_d_px.plane_label] = ...
    GraphCut('expand',gch_crf_d_px);
var_assig_crf_d_px.class = var_assig_crf_d_px.class + 1;
gch_crf_d_px = GraphCut('close', gch_crf_d_px);
acc_crf_d_px = sum(var_assig_crf_d_px.class == Label) / numel(Label);

% Crf+D+PX+TP model
disp('==========================================');
disp('Starting Crf+D+PX+TP model');
disp('==========================================');
gch_crf_d_px_tp = GraphCut('open_hop', Feat.StuffUnary, Feat.Sc, Feat.Ss, struct([]), ...
	Feat.HopDet, struct([]), struct([]), Feat.HopPlane, Feat.ThingPlane);
[gch_crf_d_px_tp var_assig_crf_d_px_tp.class ...
    var_assig_crf_d_px_tp.det_label var_assig_crf_d_px_tp.plane_label] = ...
    GraphCut('expand',gch_crf_d_px_tp);
var_assig_crf_d_px_tp.class = var_assig_crf_d_px_tp.class + 1;
gch_crf_d_px_tp = GraphCut('close', gch_crf_d_px_tp);
acc_crf_d_px_tp = sum(var_assig_crf_d_px_tp.class == Label) / numel(Label);


