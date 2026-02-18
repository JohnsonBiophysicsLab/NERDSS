for i in rev_2D rev_3D rev_3Dto2D
do
    mkdir data_$i
    cp $i/*.mol $i/trajectory.xyz $i/system.psf $i/*.inp $i/*.pdb data_$i
done