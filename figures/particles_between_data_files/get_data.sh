#!/bin/bash

rsync kawak@circe.rc.usf.edu:/work_bgfs/k/kawak/data/extend_filled_95xlink/s3_deform_longer/traj_cutup/{all,bw,filler,polymer}*.stats .
rsync kawak@circe.rc.usf.edu:/work_bgfs/k/kawak/data/extend_filled_95xlink/s3_deform_longer/dataframe.out .

