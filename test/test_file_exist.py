import os
import subprocess


run = 'test_file_exist'

if run == 'test_file_exist':
    prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects'
    subjects = os.listdir(prepath)
    pospath = "ROI/benson_retinotopic_map"

    pathtck = [os.path.join(prepath,subject,pospath) for subject in subjects]
    tck_exist = map(os.path.exists,pathtck)

    cd = dict(zip(subjects,tck_exist))
    for k in sorted(cd):
        print("{}:{}".format(k,cd[k]))

elif run == 'subjects_list':
    # list subjects in the special path
    prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects'
    subjects = os.listdir(prepath)
    subjects = sorted(subjects)
    subjects.pop()
    print(subjects)
    subjectss = map(int,subjects)
    print(subjectss)
    for i in subjectss:
        print(i)

elif run == 'create directory batch':
# create directory in many subjects directory
    prepath = '/nfs/s2/userhome/quyukun/workingdir/fiberdata/subjects'
    subjects = os.listdir(prepath)
    subjects = sorted(subjects)
    subjects.pop()
    subjects = ['101915', '102008',
                '102109', '102311', '102513', '102614', '102715', '102816', '103010', '103111', '111211', '200008',
                '201818',  '530635', '555651', '597869', '686969', '767464', '888678', '995174', '996782']
    for subject in subjects:
        directory = os.path.join(prepath,subject,"Diffusion/tractography/endpoints")
        os.mkdir(directory)
