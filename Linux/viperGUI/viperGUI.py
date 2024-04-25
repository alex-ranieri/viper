import os
import sys
import subprocess
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QPushButton, QFileDialog, QComboBox, QMessageBox, QVBoxLayout, QWidget
from PyQt5.QtGui import QPixmap, QColor, QPalette
from PyQt5.QtCore import Qt

class ViperGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("VIPER Runner")

        # Acrylic layout
        self.setWindowOpacity(0.95)  # Adjust alpha value for transparency

        # Expand the user's home directory in the logo path
        self.logo_path = os.path.expanduser("~/viperGUI/VIPER_ico.png")

        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        self.folder_path = QLineEdit(self)
        self.num_samples = QLineEdit(self)
        self.pipeline_choice = QComboBox(self)
        self.num_threads = QComboBox(self)
        self.num_threads.setEditable(True)  # Make the combo box editable
        self.num_threads.setEditText("1")  # Default number of threads is set to 1

        # Create and set the logo
        self.logo = QLabel(self)
        self.logo.setAlignment(Qt.AlignTop)

        # Create and set up widgets
        self.create_widgets()

    def create_widgets(self):
        layout = QVBoxLayout(self.centralWidget())

        # Logo
        self.update_logo_size()  # Update logo size initially
        layout.addWidget(self.logo)

        # Folder Selection
        label_folder = QLabel("Select Folder:", self)
        layout.addWidget(label_folder)
        layout.addWidget(self.folder_path)
        button_browse = QPushButton("Browse", self)
        button_browse.clicked.connect(self.browse_folder)
        layout.addWidget(button_browse)

        # Number of Samples
        label_samples = QLabel("Number of Samples to Analyze in Parallel:", self)
        layout.addWidget(label_samples)
        layout.addWidget(self.num_samples)

        # Viral Pipeline Selection
        label_pipeline = QLabel("Select Viral Pipeline:", self)
        layout.addWidget(label_pipeline)
        layout.addWidget(self.pipeline_choice)
        self.pipeline_choice.addItems(["VIPER_CoV.sh", "VIPER_DENV.sh", "VIPER_Influenza.sh"])
        self.pipeline_choice.currentIndexChanged.connect(self.on_pipeline_changed)

        # Number of Threads Selection
        label_threads = QLabel("Number of Threads per Sample:", self)
        layout.addWidget(label_threads)
        layout.addWidget(self.num_threads)
        available_threads = os.cpu_count() or 1
        thread_values = [str(i) for i in range(1, available_threads + 1)]
        self.num_threads.addItems(thread_values)
        self.num_threads.currentIndexChanged.connect(self.on_threads_changed)

        # Run Button
        button_run = QPushButton("Run Pipeline", self)
        button_run.clicked.connect(self.run_pipeline)
        layout.addWidget(button_run)

        # Update Virus Databases Button
        button_update_databases = QPushButton("Update Virus Databases", self)
        button_update_databases.clicked.connect(self.update_databases)
        layout.addWidget(button_update_databases)

        # Apply style to adjust font size
        self.setStyleSheet(
            """
            QLabel {
                font-size: 14px;
            }
            QPushButton {
                font-size: 14px;
                padding: 8px;
            }
            """
        )

        # Apply Fusion style
        app.setStyle("Fusion")
        dark_palette = QPalette()
        dark_palette.setColor(QPalette.Window, QColor(41, 41, 41))
        dark_palette.setColor(QPalette.WindowText, Qt.white)
        dark_palette.setColor(QPalette.Base, QColor(25, 25, 25))
        dark_palette.setColor(QPalette.AlternateBase, QColor(41, 41, 41))
        dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
        dark_palette.setColor(QPalette.ToolTipText, Qt.white)
        dark_palette.setColor(QPalette.Text, Qt.white)
        dark_palette.setColor(QPalette.Button, QColor(41, 41, 41))
        dark_palette.setColor(QPalette.ButtonText, Qt.white)
        dark_palette.setColor(QPalette.BrightText, Qt.red)
        dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.HighlightedText, Qt.black)
        app.setPalette(dark_palette)

    def browse_folder(self):
        folder_selected = QFileDialog.getExistingDirectory(self, "Select Folder")
        self.folder_path.setText(folder_selected)

    def run_pipeline(self):
        folder = self.folder_path.text()
        num_samples = self.num_samples.text()
        pipeline = self.pipeline_choice.currentText()
        num_threads = self.num_threads.currentText()

        if not folder or not num_samples or not pipeline:
            QMessageBox.critical(self, "Error", "Please fill in all fields.")
            return

        try:
            # Change the working directory to the selected folder
            os.chdir(folder)
            # Run the pipeline script with the specified number of samples and threads
            command = f"{pipeline} {num_threads} {num_samples}"
            print(f"Running command: {command}")
            subprocess.run(command, check=True, shell=True)

            QMessageBox.information(self, "Executed", "VIPER executed!")
        except subprocess.CalledProcessError:
            QMessageBox.critical(self, "Error", "An error occurred while trying to run VIPER.")

    def update_databases(self):
        # Run the script to update virus databases
        pipeline_path = os.environ.get("PIPELINE")
        script_path = os.path.join(pipeline_path, "update_database", "update_virusDB.sh")
        try:
            subprocess.run([script_path], check=True, shell=True)
            QMessageBox.information(self, "Done", "Virus databases updated successfully!")
        except subprocess.CalledProcessError:
            QMessageBox.critical(self, "Error", "An error occurred while updating virus databases.")

    def on_pipeline_changed(self, index):
        selected_pipeline = self.pipeline_choice.itemText(index)
        print(f"Selected Pipeline: {selected_pipeline}")

    def on_threads_changed(self, index):
        selected_threads = self.num_threads.itemText(index)
        print(f"Selected Threads: {selected_threads}")

    def resizeEvent(self, event):
        # Called when the window is resized
        self.update_logo_size()

    def update_logo_size(self):
        # Update logo size based on window width
        logo_width = min(self.width() - 20, 500)  # Limit maximum width to 500 pixels
        self.logo.setPixmap(QPixmap(self.logo_path).scaledToWidth(logo_width, Qt.SmoothTransformation))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ViperGUI()
    window.setGeometry(100, 100, 520, 230)  # Set window size
    window.show()
    sys.exit(app.exec_())
