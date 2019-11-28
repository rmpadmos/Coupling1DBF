# 1D Blood Flow

The 1-D blood flow model generates a patient from patient-specific data, literature values and high-resolution data from volunteers.

# Blood Flow Solver
Compile the BF Solver in the Bloodflow_model folder by running msbuild or using monodevelop.
The executable location is located in Bloodflow_model/Blood Flow Model/bin/Debug/
(If necessary, run chmod +x BloodflowModel.exe)

# Patient Generation
In the python folder, run GenerateBloodflowFiles.py with executable location and patient folder location.
Example: ./GenerateBloodflowFiles.py ../Bloodflow_model/Blood Flow Model/bin/Debug/ ../Generated_Patients/patient_0/

# Running the 1-D Blood flow solver
Run BloodflowOnlyClot.py with executable location, patient folder location, and an boolean argument whether the clot should be included.
Note that the script folder location is needed by the bf solver and is listed in Model_parameters.txt. The location is updated during the patient generation but can be different if the patient is generated on another machine.
Example: ./BloodflowOnlyClot.py ../Bloodflow_model/Blood Flow Model/bin/Debug/BloodflowModel.exe ../Generated_Patients/patient_0/ "True"

