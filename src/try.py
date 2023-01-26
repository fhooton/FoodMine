"""import os

os.system("jupyter notebook")

docker run  -v  D:\m:/MyProjectFiles/in  -v D:\o:/MyProjectFiles/out brain-docker

os.system("some_command < input_file | another_command > output_file")

os.system("docker run  -v < D:\m:/MyProjectFiles/in | -v D:\o:/MyProjectFiles/out > D:\o:/MyProjectFiles/out ")
"""

import subprocess

with open("/tmp/output.log", "a") as output:
    subprocess.call("docker run  -v  D:\m:/MyProjectFiles/in  -v D:\o:/MyProjectFiles/out brain-docker", shell=True,
                    stdout=output, stderr=output)
