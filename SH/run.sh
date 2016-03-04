# step 1, make vs and density model
#matlab -nodesktop -nosplash -nodisplay <<END
#make_prem_model_sh
#exit
#END

# step 2, calculate first arrival travel time
#punch par=prem_withpunch.par



# step 3, run finite difference
#nbsh2d par=prem_withpunch.par
#
## step 4, get sacfile
matlab -nodesktop -nosplash -nodisplay <<END
nt=120001;
isis2vel_true_pos('slab120_V.isis',nt);
exit
END
