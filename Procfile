#################### How to run HostPhinder command line ######################
# From folder 
# /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder
./scripts/wrap_hostPhinder_cl.sh -i <input phage genome> -o <output directory> \
-t: <taxonomy [species, genus]>

#################### Run HostPhinder on queueing system ######################
xmsub -W group_list=pr_phage -A pr_phage \
-l nodes=1:ppn=2,mem=10gb,walltime=10:00 \
-d /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder \
-N D9 \ # Change to desired job name
-de \
-e /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder/D9.err \ # Change to desired error filename
-o /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder/D9.out \ # Change to desired output filename
-r y scripts/wrap_hostPhinder_clNEW.sh \
-i /home/projects/pr_phage/people/mvla/georgian_phages/data/D9_allContigs.fasta -t species \
-o /home/projects/pr_phage/people/juliav/HostPhinder_general/git_HP_repo/HostPhinder # change to desired output folder

