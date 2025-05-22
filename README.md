# Jupiter Lagrange Point Solar System Simulation

This project is a POC of measuring time dilation between the 5 Lagrange points of Jupiter as it  
orbits around the Sun.

## Prerequisites

* **Python 3.8+** installed on your system. You can download it from [python.org](https://www.python.org/downloads/).  
* **Git** for cloning the repository.  
* **pytest** is included in `requirements.txt` for running the test suite.

## Setup Instructions

Follow these steps to set up a virtual environment and install project dependencies:

1. **Clone the repository**  

```bash
   git clone <repo_url>
   cd <repo_dir>
````

2. **Create a virtual environment**

```bash
python3 -m venv .venv
```

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

With the virtual environment activated, run:

```bash
python src/main.py
```

Replace `main.py` with the entry point of your project.

## Running the Test Suite

This project uses **pytest**. Once your environment is set up:

1. **Ensure pytest is installed** (should be in `requirements.txt`):

   ```bash
   pip install pytest
   ```

2. **Run all tests** from the project root:

   ```bash
   pytest
   ```

   Or to run a specific test file:

   ```bash
   pytest tests/test_planet.py
   pytest tests/test_satellite.py
   ```

3. **For more verbose output**:

   ```bash
   pytest -v
   ```

## Deactivating the Virtual Environment

When you're done working:

```bash
deactivate
```

## .gitignore

Ensure you have the virtual environment folder listed in your `.gitignore`:

```gitignore
# Virtual environment
.venv/
```

## License

This project is released under the [MIT License](LICENSE).
See the `LICENSE` file for full details.

---