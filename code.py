"""
The construction of initial dataset used in DPGEN itertion by the way of CI-NEB combined with pertubation.

Author: Taiping Hu
Email: taipinghu@pku.edu.cn
"""

import json
from typing import List
from dflow import (
    Workflow,
    Step,
    argo_range,
    SlurmRemoteExecutor,
    upload_artifact,
    download_artifact,
    InputArtifact,
    OutputArtifact,
    ShellOPTemplate
)
from dflow.python import (
    PythonOPTemplate,
    OP,
    OPIO,
    OPIOSign,
    Artifact,
    Slices
)

import subprocess, os, shutil, glob
from pathlib import Path
from typing import List
from dflow import config, s3_config


def create_path(path):
    path += '/'
    if os.path.isdir(path) : 
        dirname = os.path.dirname(path)        
        counter = 0
        while True :
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname) : 
                shutil.move (dirname, bk_dirname) 
                break
                counter += 1
    os.makedirs (path)
    
def select_logs(root_path, fname='OUTCAR'):
    logs = []
    for r, d, fs in os.walk(root_path):
        for f in fs:
            if os.path.islink(f):
                continue
            if f == fname:
                logs.append(os.path.join(r, f))
    return sorted(logs)

class VASPOpt(OP):
    """
    class for VASP optimization
    """
    def __init__(self, 
                intel_env: str,   # intel environment path
                vasp_exec_path: str,  # vasp execute path
                ncores: int # number of cores used in the calculations
                ):
        self.__intel_env = intel_env
        self.__vasp_exec_path = vasp_exec_path
        self.__ncores = ncores

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            'input': Artifact(Path),
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "contcar": Artifact(Path),
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        cwd = os.getcwd()
        os.chdir(op_in["input"])
        cmd = "ulimit -s unlimited; source %s; mpiexec.hydra -genv I_MPI_DEBUG 0 -genv I_MPI_DEVICE ssm -np %d %s"%(self.__intel_env, self.__ncores, self.__vasp_exec_path)
        subprocess.call(cmd, shell=True)
        os.chdir(cwd)
        op_out = OPIO({
            "contcar": Path(op_in["input"])/"CONTCAR",  # return CONTCAR for next CI-NEB step
        })
        return op_out

class CINEB(OP):
    """
    perform CI-NEB to construct intermediate images
    """
    def __init__(self, vtst_script_path):
        self.__vtst_script_path = vtst_script_path

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "image_numb": int,
            "contcar": Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "neb_path": Artifact(Path)
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        cwd = os.getcwd()
        os.chdir(op_in["contcar"])
        cmd = self.__vtst_script_path + "/nebmake.pl initial/CONTCAR final/CONTCAR %d"%op_in["image_numb"]
        subprocess.call(cmd, shell=True)
        os.system("mkdir CINEB; mv 0* CINEB")
        os.chdir(cwd)

        op_out = OPIO({
            "neb_path": op_in["contcar"]/"CINEB"  
        })
        return op_out

class Pertubation(OP):
    """
    Pertubation for intermediate images 
    """
    def __init__(self, pert_script_path):
        self.__pert_script_path = pert_script_path

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "pert_numb": int,
            "incar_cineb": Artifact(Path),
            "kpoints": Artifact(Path),
            "potcar": Artifact(Path),
            "neb_path": Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "pert_path": Artifact(List[Path])
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        neb_numb = 3
        cwd = os.getcwd()
        neb_path = os.path.abspath(str(op_in["neb_path"]))
        
        INCAR = os.path.abspath(str(op_in["incar_cineb"]))
        KPOINTS = os.path.abspath(str(op_in["kpoints"]))
        POTCAR = os.path.abspath(str(op_in["potcar"]))

        final_dirs = []
        for ii in range(op_in["pert_numb"]):
            os.chdir(op_in["neb_path"])
            tmp_path = os.path.abspath(os.path.join("./", "Pertubation/Pert-%02d"%ii))
            final_dirs.append(Path(tmp_path))
            create_path(tmp_path)
            os.system("cp %s Pertubation/Pert-%02d/INCAR"%(INCAR, ii))
            os.system("cp %s Pertubation/Pert-%02d/POTCAR"%(POTCAR, ii))
            os.system("cp %s Pertubation/Pert-%02d/KPOINTS"%(KPOINTS, ii))
            os.chdir(tmp_path)
            for jj in range(neb_numb+2):
                neb_dir = "%02d"%jj
                poscar_dir = os.path.join(neb_path, neb_dir)
                poscar = os.path.join(poscar_dir, 'POSCAR')
                work_path = os.path.join(tmp_path, neb_dir)
                create_path(work_path)
                if jj == 0 or jj == neb_numb+1:
                    pos_in = poscar
                    pos_out = os.path.join(work_path, 'POSCAR')
                    os.system("cp %s %s "%(pos_in, pos_out))
                else:
                    pert_cmd = "python3 " + self.__pert_script_path + ' -etmax %f -ofmt vasp %s %d %f > /dev/null' %(1, poscar, 1, 0.01)
                    subprocess.check_call(pert_cmd, shell=True)
                    pos_in = os.path.join(poscar_dir, 'POSCAR1.vasp')
                    pos_out = os.path.join(work_path, 'POSCAR')
                    os.system("cp %s %s "%(pos_in, pos_out))
                    os.remove(pos_in)
            os.chdir(cwd)

        op_out = OPIO({
            "pert_path": final_dirs
        })
        return op_out

class VASPVTSTOpt(OP):
    """
    Optimization using VASP VTST version
    """
    def __init__(self, 
                intel_env: str,   # intel environment path
                vasp_exec_path: str,  # vasp execute path
                ncores: int # number of cores used in the calculations
                ):
        self.__intel_env = intel_env
        self.__vasp_exec_path = vasp_exec_path
        self.__ncores = ncores

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "pert_path": Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "final_path": Artifact(Path)
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        print(" Current work dir is %s"%op_in["pert_path"])
        os.chdir(op_in["pert_path"])
        cmd = "ulimit -s unlimited; source %s; mpiexec.hydra -genv I_MPI_DEBUG 0 -genv I_MPI_DEVICE ssm -np %d %s"%(self.__intel_env, self.__ncores, self.__vasp_exec_path)
        subprocess.call(cmd, shell=True)
        op_out = OPIO({
            "final_path": op_in["pert_path"]
        })
        return op_out

class DataCollection(OP):
    """
    Data collection using dpdata
    """
    def __init__(self):
        pass

    @classmethod
    def get_input_sign(cls):
        return OPIOSign({
            "final_path": Artifact(Path)
        })

    @classmethod
    def get_output_sign(cls):
        return OPIOSign({
            "data_dir": Artifact(Path)
        })

    @OP.exec_sign_check
    def execute(self, op_in: OPIO) -> OPIO:
        import dpdata

        final_sys = None

        cwd = os.getcwd()
        os.chdir(op_in["final_path"])
        all_outcars = select_logs("./", "OUTCAR")

        for outcar in all_outcars:
            ls = dpdata.LabeledSystem(outcar, fmt="vasp/outcar")
            for iframe in range(ls.get_nframes()):
                if final_sys is None:
                    final_sys = ls[iframe]
                else:
                    final_sys.append(ls[iframe])
        final_sys.to_deepmd_raw("deepmd")
        final_sys.to_deepmd_npy("deepmd")

        os.chdir(cwd)

        op_out = OPIO({
            "data_dir": op_in["final_path"] / "deepmd"
        })
        return op_out

def main():
    slurm_remote_executor = SlurmRemoteExecutor(host="", port=22, username="",  password="", header="#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=128\n#SBATCH --job-name=VASP\n#SBATCH --partition=amd_512\n", pvc=None)

    jdata = json.load(open("pert.json"))

    Vasp_Opt = Step(
        "VASP-Opt",
        PythonOPTemplate(VASPOpt(jdata["intel_env_path"], jdata["vasp_exec_path"], jdata["ncores"]), image="dptechnology/dflow", slices=Slices("{{item}}", input_artifact=["input"], output_artifact=["contcar"])),
        artifacts={"input": upload_artifact([jdata["initial_state"], jdata["final_state"]])},
        with_param=range(2),
        key="VASPOpt-{{item}}",
        executor=slurm_remote_executor,
    )

    CI_NEB = Step(
        "CI-NEB-Generation",
        PythonOPTemplate(CINEB(jdata["vtst_script_path"], ), image="dptechnology/dflow"),
        parameters={"image_numb": jdata["image_numb"]},
        artifacts={"contcar": Vasp_Opt.outputs.artifacts["contcar"]},
        executor=slurm_remote_executor,
    )

    Pert = Step(
        "Pertubation",
        PythonOPTemplate(Pertubation(jdata["pert_script_path"]), image="dptechnology/dflow"),
        parameters={"pert_numb": jdata["pert_numb"]},
        artifacts={"neb_path": CI_NEB.outputs.artifacts["neb_path"], "incar_cineb":upload_artifact([jdata["incar_neb"]]), "kpoints":upload_artifact([jdata["kpoints"]]), "potcar":upload_artifact([jdata["potcar"]])},
        executor=slurm_remote_executor,
    )

    VASP_VTST_Opt = Step(
        "VASP-VTST-Opt",
        PythonOPTemplate(VASPVTSTOpt(jdata["intel_env_path"], jdata["vasp_exec_path"], jdata["ncores"]), image="dptechnology/dflow", slices=Slices("{{item}}", input_artifact=["pert_path"], output_artifact=["final_path"])),
        artifacts={"pert_path": Pert.outputs.artifacts["pert_path"]},
        with_param=range(jdata["pert_numb"]),
        key="VASPVTSTOpt-{{item}}",
        executor=slurm_remote_executor,
    )

    Data_Collection = Step(
        "Data-Collection",
        PythonOPTemplate(DataCollection, image="dptechnology/dflow"),
        artifacts={"final_path": VASP_VTST_Opt.outputs.artifacts["final_path"]},
        executor=slurm_remote_executor,
    )

    wf = Workflow("vasp-cineb-pert-initial-data")
    wf.add(Vasp_Opt)
    wf.add(CI_NEB)
    wf.add(Pert)
    wf.add(VASP_VTST_Opt)
    wf.add(Data_Collection)
    wf.submit()

if __name__ == "__main__":
    main()

