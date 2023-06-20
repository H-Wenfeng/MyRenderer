<div align="center">
  <a href="https://github.com/H-Wenfeng/MyRenderer">
    <img src="./2023-06-20 22-01-42[00_00_01--00_00_19].gif" width="200" height="200" />
  </a>

  <h1>ΩRenderer</h1>

  <p>
 Aiming to build a pure self-made C++-based renderer.
  </p>
<div align="left">

<!-- ![img](https://github.com/H-Wenfeng/MyRenderer/blob/main/2023-06-19%2012-21-39%5B00_00_00--00_00_20%5D.gif) -->

## Introductions

This project aims to document my self-learning journey in Computer Graphics and improve it further in the future. One of the main goals is to ensure that all important code blocks are accompanied by detailed annotations for better understanding. Additionally, I will strive to add more functions and reproduce interesting papers using this renderer.
## Features
<div align="left">
<table rules="none">
<tr>
<td>
<p>* Shader Based</p>
<p>* Depth Buffer</p>
<p>* Interpolation Normal Vector</p>
<p>* Normal Mapping  </p>
<p>* Tangent Space Normal Mapping</p>
<p>* Shadow Mapping</p>
<p>* Blinn-Phong Model</p>
<p>* X11 API + Keyboard-based Camera/Light control<p>
</td>
<td>
<right> <img src="./2023-06-19 12-21-39[00_00_00--00_00_20].gif" width="200" height="200" /></td>
</table>
</td>
</tr>
</div>

## Usage
This project is powered by WSL2 with Ubuntu 20.04.

```
sudo apt install gcc libx11-dev
make ./main
Use WASDIK to control the camera.
Use 8246+- to control the light.
```
## Reference
https://nostarch.com/computer-graphics-scratch

https://learnopengl.com/

https://github.com/ssloy/tinyrenderer
