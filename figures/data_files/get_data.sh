#!/bin/bash

prefix=filled_deform_
rsync kawak@circe.rc.usf.edu:/work_bgfs/k/kawak/data/extend_filled_95xlink/s3_deform_longer/dataframe.* .
mv dataframe.out ${prefix}dataframe.out; mv dataframe.binned.out ${prefix}dataframe.binned.out 

prefix=filled_npt_
rsync kawak@circe.rc.usf.edu:/work_bgfs/k/kawak/data/extend_filled_95xlink/s3_npt_samestyle_as_deform_longer/dataframe.* .
mv dataframe.out ${prefix}dataframe.out; mv dataframe.binned.out ${prefix}dataframe.binned.out 

prefix=neat_deform_
rsync kawak@circe.rc.usf.edu:/work_bgfs/k/kawak/data/extend_neat_95xlink/s3_deform_longer/dataframe.* .
mv dataframe.out ${prefix}dataframe.out; mv dataframe.binned.out ${prefix}dataframe.binned.out 

prefix=neat_npt_
rsync kawak@circe.rc.usf.edu:/work_bgfs/k/kawak/data/extend_neat_95xlink/s3_npt/dataframe.* .
mv dataframe.out ${prefix}dataframe.out; mv dataframe.binned.out ${prefix}dataframe.binned.out 
