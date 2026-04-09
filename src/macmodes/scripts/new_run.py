import os
import sys
import shutil
from macmodes.constants import MACMODES_DIR

def main():
    '''makes new run dir with editable config file, is named run_dir_name if passed
    as command line argument, else names it template'''

    if len(sys.argv) > 1:
       run_dir_name=sys.argv[1] 
    else:
        run_dir_name='template_run'

    template_path = os.path.join(MACMODES_DIR,"templates", "template_run")
    dest_path = os.path.join(MACMODES_DIR, "runs", run_dir_name)
    #os.mkdir(os.path.dirname(dest_path)) 
    shutil.copytree(template_path,dest_path) 

if __name__=="__main__":
    main()
