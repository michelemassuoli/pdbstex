# Troubleshooting: installing and running `pdbstex` (Windows + Linux)

This guide fixes the most common issues where `pip install pdbstex` appears to succeed but the `pdbstex` command is not found, and where Linux shows **externally-managed-environment (PEP 668)**.

> Golden rule (everywhere): prefer `python -m pip ...` over `pip ...` when working inside environments (venv/conda/mamba). That guarantees pip belongs to the Python you are using.

---

## 0) Quick diagnosis (Windows and Linux)

Run these commands and check the paths:

### Check which Python/pip you are using
**Windows (cmd/PowerShell)**
```bat
where python
where pip
python -c "import sys; print(sys.executable)"
python -m pip -V
```

**Linux (bash)**
```bash
which python
which pip
python -c "import sys; print(sys.executable)"
python -m pip -V
```

### Check the `pdbstex` command
**Windows**
```bat
where pdbstex
pdbstex -h
```

**Linux**
```bash
which pdbstex
pdbstex -h
```

If `python`/`pip` point to the system Python instead of your environment, you are in the “environment name shown, but not actually active” case (see below).

---

## 1) Windows

### 1.1 Common issue: prompt shows `(prova)` but `python` comes from Windows Store (`WindowsApps`)
Symptoms:
- `where python` shows `C:\Users\...\AppData\Local\Microsoft\WindowsApps\python.exe`
- `pip install pdbstex` installs into `...LocalCache\local-packages...`
- `pdbstex -h` prints: “not recognized as an internal or external command”

**Why it happens**
You are using the Microsoft Store shims in PATH, not the venv/conda Python. As a result, the `pdbstex` launcher is not created in your environment (or is not on PATH).

#### Fix A (recommended): activate the venv for real
If your venv is called `prova` in the current directory:

**cmd.exe**
```bat
prova\Scripts\activate.bat
where python
python -c "import sys; print(sys.executable)"
python -m pip install --upgrade pip
python -m pip install --force-reinstall pdbstex
where pdbstex
pdbstex -h
```

**PowerShell**
```powershell
.\prova\Scripts\Activate.ps1
where python
python -c "import sys; print(sys.executable)"
python -m pip install --upgrade pip
python -m pip install --force-reinstall pdbstex
where pdbstex
pdbstex -h
```

If PowerShell blocks scripts:
```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

#### Fix B (no activation needed): call the venv Python directly
```bat
prova\Scripts\python.exe -m pip install --upgrade pip
prova\Scripts\python.exe -m pip install --force-reinstall pdbstex
prova\Scripts\pdbstex.exe -h
```

#### Fix C (if you intentionally used `--user` installs): add your user Scripts directory to PATH
If pip installs into user site, launchers land in a `...\Scripts` folder.
Find it with:
```bat
python -m site --user-base
```
Then add the `Scripts` directory under that base to your user PATH.

---

### 1.2 Installing inside Conda (Anaconda/Miniconda/Mambaforge)
1) Create an environment and install:
```bat
conda create -n pdbstex-env python=3.13 -y
conda activate pdbstex-env
python -m pip install --upgrade pip
python -m pip install pdbstex
pdbstex -h
```

2) Confirm you are using the env Python:
```bat
where python
python -m pip -V
```

---

## 2) Linux (Ubuntu/Debian): **externally-managed-environment (PEP 668)**

### 2.1 What it means
You are trying to install packages with `pip` into the system Python managed by `apt`. Ubuntu/Debian block this to prevent conflicts.

**Correct solutions**: use `venv`, `pipx`, or a conda/mamba environment.

---

### 2.2 Solution A: `venv` (best for a project/library workflow)
If `python -m venv` fails because `ensurepip` is missing, install the venv package:

**For Python 3.13 specifically**
```bash
sudo apt update
sudo apt install -y python3.13-venv
```

**Or generic**
```bash
sudo apt update
sudo apt install -y python3-venv
```

Then create and use the venv:
```bash
rm -rf prova
python3 -m venv prova
source prova/bin/activate
python -m pip install --upgrade pip
python -m pip install pdbstex
pdbstex -h
```

---

### 2.3 Solution B: `pipx` (best for a global CLI install)
```bash
sudo apt update
sudo apt install -y pipx
pipx ensurepath
```
Restart your shell (or `source ~/.profile`), then:
```bash
pipx install pdbstex
pdbstex -h
```

Upgrade later:
```bash
pipx upgrade pdbstex
```

---

## 3) Linux: Conda / Micromamba (and why it can “look active” but still use system Python)

### 3.1 Common symptom: `micromamba activate prova` but `python` is still `/usr/bin/python`
If inside `(prova)` you get:
```bash
which python
# -> /usr/bin/python
python -m pip -V
# -> .../dist-packages/pip ...
```
then you are **not** using the environment’s Python.

Common reasons:
- the env does not contain `python`/`pip`
- micromamba is not initialized for your shell, so activation does not prepend the env `bin/` to PATH

---

### 3.2 Robust fix (recommended): use `micromamba run`
1) Ensure Python and pip exist inside the env:
```bash
micromamba install -n prova -c conda-forge python pip
```

2) Install and run using the env explicitly:
```bash
micromamba run -n prova python -m pip install --upgrade pip
micromamba run -n prova python -m pip install pdbstex
micromamba run -n prova pdbstex -h
```

3) Verify the Python path:
```bash
micromamba run -n prova which python
micromamba run -n prova python -c "import sys; print(sys.executable)"
```

---

### 3.3 Fix activation: initialize micromamba for your shell (bash)
```bash
micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
micromamba activate prova
which python
```

After correct activation, `which python` should point to something like:
`.../micromamba/envs/prova/bin/python`

---

### 3.4 Create a clean env (recommended if `prova` is messy)
```bash
micromamba create -n pdbstex-env -c conda-forge python=3.13 pip -y
micromamba run -n pdbstex-env python -m pip install --upgrade pip
micromamba run -n pdbstex-env python -m pip install pdbstex
micromamba run -n pdbstex-env pdbstex -h
```

---

## 4) `pdbstex`-specific checks

### 4.1 Check installed version
```bash
python -c "import pdbstex; print(pdbstex.__version__)"
```

### 4.2 Show installation location
```bash
python -m pip show pdbstex
```

### 4.3 Check whether the entry point is on PATH
**Linux**
```bash
python -c "import shutil; print(shutil.which('pdbstex'))"
```

**Windows**
```bat
python -c "import shutil; print(shutil.which('pdbstex'))"
```

---

## 5) Final checklist (if `pdbstex` still won’t run)
1) `python -m pip -V` must show the environment Python (venv/conda/mamba), not system Python.
2) Reinstall inside the environment:
   - `python -m pip install --force-reinstall pdbstex`
3) Confirm the launcher exists:
   - `where pdbstex` (Windows) / `which pdbstex` (Linux)
4) If you see PEP 668 on Ubuntu/Debian:
   - use `venv`, `pipx`, or `micromamba run` (not system pip).
