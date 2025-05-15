# Jupiter Lagrange Point Solar System Simulation

This project is a POC of measuring time dilation between the 5 lagrange points of Jupiter as it
orbits around the sun.


## Prerequisites

* **Python 3.8+** installed on your system. You can download it from [python.org](https://www.python.org/downloads/).
* **Git** for cloning the repository.

## Setup Instructions

Follow these steps to set up a virtual environment and install project dependencies:

1. **Clone the repository**

   ```bash
   git clone <repo_url>
   cd <repo_dir>
   ```

2. **Create a virtual environment**

   ```bash
   python3 -m venv .venv
   ```

   This will create a folder named `.venv/` in your project directory.

3. **Activate the virtual environment**

   * On **macOS/Linux**:

     ```bash
     source .venv/bin/activate
     ```

   * On **Windows (PowerShell)**:

     ```powershell
     .\.venv\Scripts\Activate.ps1
     ```

   * On **Windows (cmd.exe)**:

     ```cmd
     .\.venv\Scripts\activate.bat
     ```

4. **Upgrade pip** (optional but recommended)

   ```bash
   pip install --upgrade pip
   ```

5. **Install project dependencies**

   ```bash
   pip install -r requirements.txt
   ```

6. **Verify the installation**

   ```bash
   python --version
   pip freeze
   ```

   You should see your Python version and a list of installed packages matching those in `requirements.txt`.

## Running the Application

With the virtual environment activated, you can run your application:

```bash
python main.py
```

Replace `main.py` with the entry-point of your project.

## Deactivating the Virtual Environment

When you're done working, deactivate the environment:

```bash
deactivate
```

## .gitignore

Ensure you have the virtual environment folder listed in your `.gitignore` to avoid committing it:

```
# Virtual environment
.venv/
```

---

Happy coding!
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
