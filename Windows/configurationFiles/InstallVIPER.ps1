 # Function to create GUI for selecting installation folder and tar.gz file
    function Show-WSL-Import-GUI {
        Add-Type -AssemblyName System.Windows.Forms
        Add-Type -AssemblyName System.Drawing

        # Create form
        $form = New-Object Windows.Forms.Form
        $form.Text = "WSL 2 VIPER Installer"
        $form.Size = New-Object Drawing.Size(400, 200)
        $form.StartPosition = "CenterScreen"

        # Create labels
        $label1 = New-Object Windows.Forms.Label
        $label1.Location = New-Object Drawing.Point(10, 20)
        $label1.Size = New-Object Drawing.Size(200, 20)
        $label1.Text = "Select installation folder for VIPER:"
        $form.Controls.Add($label1)

        $label2 = New-Object Windows.Forms.Label
        $label2.Location = New-Object Drawing.Point(10, 70)
        $label2.Size = New-Object Drawing.Size(200, 20)
        $label2.Text = "Select Ubuntu tar.gz file:"
        $form.Controls.Add($label2)

        # Create textboxes
        $folderTextBox = New-Object Windows.Forms.TextBox
        $folderTextBox.Location = New-Object Drawing.Point(10, 40)
        $folderTextBox.Size = New-Object Drawing.Size(200, 20)
        $form.Controls.Add($folderTextBox)

        $fileTextBox = New-Object Windows.Forms.TextBox
        $fileTextBox.Location = New-Object Drawing.Point(10, 90)
        $fileTextBox.Size = New-Object Drawing.Size(200, 20)
        $form.Controls.Add($fileTextBox)

        # Create browse buttons
        $folderButton = New-Object Windows.Forms.Button
        $folderButton.Location = New-Object Drawing.Point(220, 40)
        $folderButton.Size = New-Object Drawing.Size(75, 23)
        $folderButton.Text = "Browse"
        $folderButton.Add_Click({
            $folderBrowser = New-Object Windows.Forms.FolderBrowserDialog
            $folderBrowser.Description = "Select installation folder for VIPER"
            $result = $folderBrowser.ShowDialog()
            if ($result -eq [Windows.Forms.DialogResult]::OK) {
                $folderTextBox.Text = $folderBrowser.SelectedPath
            }
        })
        $form.Controls.Add($folderButton)

        $fileButton = New-Object Windows.Forms.Button
        $fileButton.Location = New-Object Drawing.Point(220, 90)
        $fileButton.Size = New-Object Drawing.Size(75, 23)
        $fileButton.Text = "Browse"
        $fileButton.Add_Click({
            $fileBrowser = New-Object Windows.Forms.OpenFileDialog
            $fileBrowser.Filter = "tar.gz files (*.tar.gz)|*.tar.gz"
            $result = $fileBrowser.ShowDialog()
            if ($result -eq [Windows.Forms.DialogResult]::OK) {
                $fileTextBox.Text = $fileBrowser.FileName
            }
        })
        $form.Controls.Add($fileButton)

        # Create import button
        $importButton = New-Object Windows.Forms.Button
        $importButton.Location = New-Object Drawing.Point(10, 130)
        $importButton.Size = New-Object Drawing.Size(75, 25)
        $importButton.Text = "Install"
        $importButton.Add_Click({
            # Get installation folder and tar.gz file paths
            $installFolder = $folderTextBox.Text
            $tarGzFile = $fileTextBox.Text

            # Validate paths
            if (-not (Test-Path -Path $installFolder -PathType Container)) {
                $invalidFolderMessage = New-Object -ComObject wscript.shell
                $invalidFolderMessage.Popup("Invalid installation folder. Please select a valid folder.", 0, "Error", 0 + 48)
                return
            }

            if (-not (Test-Path -Path $tarGzFile -PathType Leaf)) {
                $invalidFileMessage = New-Object -ComObject wscript.shell
                $invalidFileMessage.Popup("Invalid tar.gz file. Please select a valid file.", 0, "Error", 0 + 48)
                return
            }

            # Import the WSL 2 distribution
            wsl.exe --import VIPER $installFolder $tarGzFile --version 2

            $importCompleteMessage = New-Object -ComObject wscript.shell
            $importCompleteMessage.Popup("WSL 2 VIPER import completed.", 0, "WSL 2 VIPER Import Status", 1 + 48)
            $form.Close()
        })
        $form.Controls.Add($importButton)

        # Show the form
        $form.ShowDialog()
    }

    # Run GUI function for WSL import
    Show-WSL-Import-GUI

Add-Type -AssemblyName System.Windows.Forms
Add-Type -AssemblyName System.Management

# Get the directory of the script
$scriptDirectory = Split-Path -Parent $MyInvocation.MyCommand.Definition

# Function to get total RAM in GB
function Get-TotalRAM {
    $totalMemory = Get-WmiObject Win32_ComputerSystem | Select-Object TotalPhysicalMemory
    $totalRAMGB = [math]::Round($totalMemory.TotalPhysicalMemory / 1GB, 2)
    return $totalRAMGB
}

# Function to get the number of threads
function Get-ThreadCount {
    return (Get-WmiObject Win32_ComputerSystem).NumberOfLogicalProcessors
}

# Create the main form
$form = New-Object System.Windows.Forms.Form
$form.Text = "VIPER configuration"
$form.Size = New-Object System.Drawing.Size(400, 250)
$form.StartPosition = "CenterScreen"

# Create a label for the password input
$passwordLabel = New-Object System.Windows.Forms.Label
$passwordLabel.Text = "Password:"
$passwordLabel.Location = New-Object System.Drawing.Point(10, 20)
$passwordLabel.Size = New-Object System.Drawing.Size(90, 20)

# Create a text box for password input with the PasswordChar set to *
$passwordTextBox = New-Object System.Windows.Forms.TextBox
$passwordTextBox.Location = New-Object System.Drawing.Point(100, 20)
$passwordTextBox.Size = New-Object System.Drawing.Size(200, 20)
$passwordTextBox.PasswordChar = '*'

# Get total RAM and thread count
$totalRAM = Get-TotalRAM
$threadCount = Get-ThreadCount

# Display total RAM and thread count in the form
$systemInfoLabel = New-Object System.Windows.Forms.Label
$systemInfoLabel.Text = "Total RAM: $totalRAM GB | Threads: $threadCount"
$systemInfoLabel.Location = New-Object System.Drawing.Point(10, 50)
$systemInfoLabel.Size = New-Object System.Drawing.Size(300, 20)

# Create a label for threads input
$threadsLabel = New-Object System.Windows.Forms.Label
$threadsLabel.Text = "Number of Threads:"
$threadsLabel.Location = New-Object System.Drawing.Point(10, 80)
$threadsLabel.Size = New-Object System.Drawing.Size(150, 20)

# Create a text box for threads input
$threadsTextBox = New-Object System.Windows.Forms.TextBox
$threadsTextBox.Location = New-Object System.Drawing.Point(180, 80)
$threadsTextBox.Size = New-Object System.Drawing.Size(50, 20)

# Create a label for RAM input
$ramLabel = New-Object System.Windows.Forms.Label
$ramLabel.Text = "Amount of RAM (GB):"
$ramLabel.Location = New-Object System.Drawing.Point(10, 110)
$ramLabel.Size = New-Object System.Drawing.Size(150, 20)

# Create a text box for RAM input
$ramTextBox = New-Object System.Windows.Forms.TextBox
$ramTextBox.Location = New-Object System.Drawing.Point(180, 110)
$ramTextBox.Size = New-Object System.Drawing.Size(50, 20)

# Create a button to execute the configuration
$executeButton = New-Object System.Windows.Forms.Button
$executeButton.Text = "Start configuration"
$executeButton.Location = New-Object System.Drawing.Point(10, 140)
$executeButton.Size = New-Object System.Drawing.Size(150, 30)
$executeButton.Add_Click({
    # Read values from text boxes
    $password = $passwordTextBox.Text
    $threads = $threadsTextBox.Text
    $ram = $ramTextBox.Text

    # Write information to .wslconfig
    $wslConfigPath = "$env:USERPROFILE\.wslconfig"
    $wslConfigContent = @"
[wsl2]
memory=${ram}GB 
processors=$threads
"@
    $wslConfigContent | Set-Content -Path $wslConfigPath -Force

    $bashCommand = """sh GenomeAssembler_configuration.sh $($password)"""

    # Start a new PowerShell process to execute the Bash command
    Start-Process wsl -ArgumentList @("-d", "VIPER", "bash", "-c", $bashCommand)
})

# Create an "Exit" button
$exitButton = New-Object System.Windows.Forms.Button
$exitButton.Text = "Exit"
$exitButton.Location = New-Object System.Drawing.Point(170, 140)
$exitButton.Size = New-Object System.Drawing.Size(150, 30)
$exitButton.Add_Click({
    # Close the form to terminate the program
    wsl.exe -t VIPER
    wsl.exe --shutdown
    $form.Close()
})

# Add controls to the form
$form.Controls.Add($passwordLabel)
$form.Controls.Add($passwordTextBox)
$form.Controls.Add($systemInfoLabel)
$form.Controls.Add($threadsLabel)
$form.Controls.Add($threadsTextBox)
$form.Controls.Add($ramLabel)
$form.Controls.Add($ramTextBox)
$form.Controls.Add($executeButton)
$form.Controls.Add($exitButton)

# Show the form
$form.ShowDialog()



