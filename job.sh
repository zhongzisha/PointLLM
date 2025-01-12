source /data/zhongz2/anaconda3/bin/activate pointllm
module load CUDA/11.8
module load cuDNN/8.9.2/CUDA-11
module load gcc/11.3.0




conda create -n pointllm python=3.9 ipython

conda activate pointllm

pip install git+https://github.com/huggingface/transformers.git@cae78c46


cd ~/DeepSpeed
mkdir -p deepspeed/ops/spatial/
export NVCC_PREPEND_FLAGS="--forward-unknown-opts"
DS_BUILD_OPS=1 DS_BUILD_CUTLASS_OPS=0 DS_BUILD_RAGGED_DEVICE_OPS=0 DS_BUILD_EVOFORMER_ATTN=0 DS_BUILD_SPARSE_ATTN=0 MAX_JOBS=16 \
pip install -e . --global-option="build_ext" --global-option="-j16" \
--cache-dir /lscratch/$SLURM_JOB_ID/deepspeed_cachedir









