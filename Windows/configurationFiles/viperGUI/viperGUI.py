import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk  # Import ttk for themed widgets
import subprocess
import os

class ViperGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("VIPER Runner")

        # Set the script folder as the current working directory
        #script_folder = os.path.dirname(os.path.realpath(__file__))
        #os.chdir(script_folder)

        # Acrylic layout
        self.master.attributes('-alpha', 0.9)  # Adjust alpha value for transparency

        self.logo_path = "~/viperGUI/VIPER_ico.png"  # Logo file located in the script folder
        self.folder_path = tk.StringVar()
        self.num_samples = tk.IntVar()
        self.pipeline_choice = tk.StringVar()
        self.num_threads = tk.StringVar(value="1")  # Default number of threads is set to 1

        # Create and set the small logo
        self.logo = tk.PhotoImage(file=self.logo_path)
        self.logo = self.logo.subsample(2, 2)  # Change the subsample values to adjust the logo size

        self.master.iconphoto(False, self.logo)

        # Create and set up widgets
        self.create_widgets()

    def create_widgets(self):
        # Logo
        tk.Label(self.master, image=self.logo).grid(row=0, column=0, columnspan=3, pady=10)

        # Folder Selection
        tk.Label(self.master, text="Select Folder:").grid(row=1, column=0, sticky="e")
        tk.Entry(self.master, textvariable=self.folder_path, width=50).grid(row=1, column=1)
        tk.Button(self.master, text="Browse", command=self.browse_folder).grid(row=1, column=2)

        # Number of Samples
        tk.Label(self.master, text="Number of Samples to Analyze in Parallel:").grid(row=2, column=0, sticky="e")
        tk.Entry(self.master, textvariable=self.num_samples).grid(row=2, column=1)

        # Viral Pipeline Selection
        tk.Label(self.master, text="Select Viral Pipeline:").grid(row=3, column=0, sticky="e")
        pipelines = ["VIPER_CoV.sh", "VIPER_DENV.sh", "VIPER_Influenza.sh"]
        tk.OptionMenu(self.master, self.pipeline_choice, *pipelines).grid(row=3, column=1)

        # Number of Threads Selection
        tk.Label(self.master, text="Number of Threads per Sample:").grid(row=4, column=0, sticky="e")
        available_threads = os.cpu_count() or 1
        thread_values = [str(i) for i in range(1, available_threads + 1)]
        ttk.Combobox(self.master, textvariable=self.num_threads, values=thread_values).grid(row=4, column=1)

        # Run Button
        tk.Button(self.master, text="Run Pipeline", command=self.run_pipeline).grid(row=5, column=0, columnspan=3, pady=10)

    def browse_folder(self):
        folder_selected = filedialog.askdirectory()
        self.folder_path.set(folder_selected)

    def run_pipeline(self):
        folder = self.folder_path.get()
        num_samples = self.num_samples.get()
        pipeline = self.pipeline_choice.get()
        num_threads = self.num_threads.get()

        if not folder or not num_samples or not pipeline:
            messagebox.showerror("Error", "Please fill in all fields.")
            return

        try:
            # Change the working directory to the selected folder
            os.chdir(folder)
            # Run the pipeline script with the specified number of samples and threads
            command = f"{pipeline} {num_threads} {num_samples}"
            print(f"Running command: {command}")
            subprocess.run(command, check=True, shell=True)

            messagebox.showinfo("Success", "VIPER executed!")
        except subprocess.CalledProcessError:
            messagebox.showerror("Error", "An error occurred while running VIPER.")

if __name__ == "__main__":
    root = tk.Tk()
    app = ViperGUI(root)
    root.mainloop()
