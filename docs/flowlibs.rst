Getting the Geographical data using FlowLibs
============================================

As part of the development of Pycequeau, an additional repository provides helper libraries for computationally
intensive hydrological terrain processing. These utilities are distributed in the following GitHub repository:

https://github.com/erinconv/FlowLibs


What is FlowLibs?
-----------------

FlowLibs provides C++ implementations of several hydrological terrain-analysis algorithms. It includes tools to compute:

- Flow direction
- Flow accumulation
- Watershed delineation
- Sub-watershed delineation

For full details (algorithmic choices, inputs/outputs, and assumptions), refer
to the official FlowLibs documentation.

How to use FlowLibs
-------------------


FlowLibs is distributed as a Docker image. Therefore, using it requires
installing Docker:

https://www.docker.com/

After installing Docker, open a terminal (or Command Prompt / PowerShell on Windows) and pull the latest compiled image:

.. code-block:: bash

    docker pull erinconv/flowlibs:latest

Next, create a Python script named ``run_carving.py`` and paste the code below. This script wraps the FlowLibs carving executable and writes outputs to a local
folder next to the input DEM.

.. code-block:: python

    import os
    import argparse
    import subprocess
    import sys


    def main() -> int:
        """Python wrapper for executing the carving algorithm from FlowLibs

        Raises:
            FileNotFoundError: If the inut file is not found

        Returns:
            int: Exit code (0 for success)
        """
        parser = argparse.ArgumentParser("Python wrapper for executing the carving algorithm from FlowLibs")

        parser.add_argument("input_path",
                            type =str,
                            help = "Path to the input DEM File")

        args = parser.parse_args()

        # Assign parsed arguments to the variablesd
        input_path = args.input_path
        if not os.path.isfile(input_path):
            raise FileNotFoundError(f"Input file not found: {input_path}", file=sys.stderr)

        mount_dir = os.path.dirname(input_path)
        input_name = os.path.basename(input_path)
        image = "erinconv/flowlibs:latest"

        # Create an output directory alongside the input
        output_dir = os.path.join(mount_dir, "dephier_outputs")
        os.makedirs(output_dir, exist_ok=True)

        # Output prefix (inside container) and ocean level
        output_prefix_in_container = "/out/run1"
        ocean_level = "0"

        # Build docker command:
        # - Mount input directory read-only at /data
        # - Mount output directory read-write at /out
        # - Invoke: dephier <command> /data/<input> /out/run1 <ocean_level> [additional args]
        cmd = [
            "docker",
            "run",
            "--rm",
            "--entrypoint", "main_carve",  # Note: no .exe extension in container
            "-v",
            f"{mount_dir}:/data:ro",
            "-v",
            f"{output_dir}:/out",
            image,
            f"/data/{input_name}",
            output_prefix_in_container,
            ocean_level,
        ]

        # if command == "2":
        #     # Add 'mod' to enable saving flow directions
        cmd.append("mod")
        print(" ".join(cmd))

        try:
            completed = subprocess.run(cmd, check=True)
            if completed.returncode != 0:
                return completed.returncode
        except subprocess.CalledProcessError as e:
            return e.returncode

        return 0


    if __name__ == "__main__":
        raise SystemExit(main())

Then execute the script as follows:

.. code-block:: bash

    python run_carving.py path/to/your/DEM.tif

Outputs
~~~~~~~

This command generates the following files in the ``dephier_outputs`` directory:

- ``run1-mod-dirs.tif``: flow direction raster
- ``run1_accu.tif``: flow accumulation raster
- ``run1-mod-dem-carved.tif``: carved (breached) DEM

Integration with the Pycequeau Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The carved DEM (``run1-mod-dem-carved.tif``) can be used as an improved elevation input for subsequent preprocessing steps (e.g., flow routing and watershed
delineation) in GRASS-GIS or other GIS software.

Using FlowLibs carving allows you to bypass the depression filling step, since depressions are treated through a breaching (carving) approach.

Alternatively, FlowLibs also provides watershed delineation utilities. Refer to the FlowLibs documentation for usage of the 
``run_watersheds.py`` workflow to produce watershed masks directly.