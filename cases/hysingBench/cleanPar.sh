#!/bin/bash
# 清理当前算例的运行产物，恢复到"已配置未运行"状态（保留 0/ constant/ system/）
# 用法:  在算例目录内 ./cleanPar.sh
[ -f system/controlDict ] || { echo "错误: 请在算例目录内运行（找不到 system/controlDict）"; exit 1; }

rm -rf processor* postProcessing log.*
find . -maxdepth 1 -type d -regextype posix-extended \
     -regex './[0-9]+(\.[0-9]+)?' ! -name 0 -exec rm -rf {} +
echo "[$(basename $PWD)] 已清理 -> $(ls | tr '\n' ' ')"
