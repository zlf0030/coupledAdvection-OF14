#!/bin/bash
# 在算例目录内运行 4 进程并行仿真
# 用法:  ./runPar.sh                        前台运行
#        nohup ./runPar.sh > log.runPar 2>&1 &   后台挂机（断终端不中断）
#        tail -f log.foamRun                看求解器进度
set -e

[ -f system/controlDict ] || { echo "错误: 请在算例目录内运行（找不到 system/controlDict）"; exit 1; }

# 加载 OpenFOAM 环境
source /home/linfan/openfoam/setup.sh >/dev/null 2>&1 \
  || source /home/linfan/openfoam/OpenFOAM-14/etc/bashrc >/dev/null 2>&1

echo "[$(date +%H:%M:%S)] === $(basename $PWD) : decomposePar ==="
decomposePar -force > log.decomposePar 2>&1

echo "[$(date +%H:%M:%S)] === $(basename $PWD) : mpirun -np 4 foamRun -parallel ==="
mpirun -np 4 foamRun -parallel > log.foamRun 2>&1
echo "[$(date +%H:%M:%S)] === $(basename $PWD) : foamRun 结束 (exit=$?) ==="

echo "[$(date +%H:%M:%S)] === $(basename $PWD) : reconstructPar ==="
reconstructPar > log.reconstructPar 2>&1
echo "[$(date +%H:%M:%S)] === $(basename $PWD) : 完成 ==="
