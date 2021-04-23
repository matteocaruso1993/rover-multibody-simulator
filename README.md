# rover-multibody-simulator package
This package simulates the motion of a rover by using exact equations of motion computed simbolically
![rover_gif](https://user-images.githubusercontent.com/70321193/112718941-6d31a080-8ef6-11eb-9c1a-13c8428d5094.gif)

Download the repository and install the package and its dependencies:
```
pip install -e .
```

To launch the demo simulator, move in the demo folder and execute:
```
python demo.py
```

Once the package has been installed, you can easily import in in any python script by doing for example:
```
from rover_multibody_simulator import four_ws_rover_dynamic_simulator as sim
```

If you want to use compiled code you need to install Microsoft C++ Build Tools. You can download it at https://visualstudio.microsoft.com/it/visual-cpp-build-tools/

Once installed the required build tools, locate in the demo folder and execute:

```
python compiler.py
```

The code will generate some compiled C code in the folder demo/data/wrappers/full-model. Once it is done you can symply run the demo_fully_parametric script:

```
python demo_fully_parametric.py
```

or you can use it as starting point and edit it for your needs.
