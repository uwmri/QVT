function [] = saveQVTData(directory,area_val,diam_val,branchList,flowPerHeartCycle_val,maxVel_val,velMean_val,nframes,matrix,res,timeres,...
    VENC,segment,PI_val,RI_val,flowPulsatile_val,r, timeMIPcrossection ,MAGcrossection,segmentFull,segmentFullJS,autoFlow,vTimeFrameave,...
    Planes,bnumMeanFlow,bnumStdvFlow,StdvFromMean,pixelSpace,VplanesAllx,VplanesAlly,VplanesAllz,imageData,caseFilePath)

    data_struct = [];
    data_struct.directory = directory;
    data_struct.area_val = area_val;
    data_struct.diam_val = diam_val;
    data_struct.branchList = branchList;
    data_struct.flowPerHeartCycle_val = flowPerHeartCycle_val;
    data_struct.maxVel_val = maxVel_val;
    data_struct.velMean_val = velMean_val;
    data_struct.nframes = nframes;
    data_struct.matrix = matrix;
    data_struct.res = res;
    data_struct.timeres = timeres;
    data_struct.VENC = VENC;
    data_struct.segment = segment;
    data_struct.PI_val = PI_val;
    data_struct.RI_val = RI_val;
    data_struct.flowPulsatile_val = flowPulsatile_val;
    data_struct.r = r;
    data_struct.timeMIPcrossection = timeMIPcrossection;
    data_struct.MAGcrossection = MAGcrossection;
    data_struct.segmentFull = segmentFull;
    data_struct.segmentFullJS = segmentFullJS; %SD
    data_struct.autoFlow = autoFlow; %SD
    data_struct.vTimeFrameave = vTimeFrameave;
    data_struct.Planes = Planes;
    data_struct.bnumMeanFlow = bnumMeanFlow;
    data_struct.bnumStdvFlow = bnumStdvFlow;
    data_struct.StdvFromMean = StdvFromMean;
    data_struct.pixelSpace = pixelSpace;
    
    Vel_Time_Res.VplanesAllx = VplanesAllx; %TR vel planes (uninterped)
    Vel_Time_Res.VplanesAlly = VplanesAlly;
    Vel_Time_Res.VplanesAllz = VplanesAllz;

    save(caseFilePath,'data_struct','Vel_Time_Res','imageData')