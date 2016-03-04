# make model, generate homo.vp homo.vs homo.den
matlab -nodesktop -nosplash <<END
mkmodel
END


# run FD
nbpsv2d par=homo.par

# convert line source to point source, generate SAC file
matlab -nodesktop -nosplash <<END
nt=2000;
isis2vel_true_pos('output_W.isis',nt);
isis2vel_true_pos('output_U.isis',nt);
exit;
END
