In the spirit of reproducible research, we provide here the MATLAB codes for replicate the results in our published article: 
**"A well-conditioned direct PinT algorithm for first- and second-order evolutionary equations"** by Jun Liu, Xiang-Sheng Wang, Shu-Lin Wu & Tao Zhou.  

*We point out the parallel results in the paper are obtained from C/MPI/Petsc implementation of the following MATLAB version algorithm.*

# Below is the short description of each code:  
  -**Compare_eigB_VDV_Table1.m**, reproduce results in Table 1  
  -**Heat1D_pint_direct_fast.m**, 1D heat equation (linear)  
  -**Heat1D_pint_direct_nonlinear.m**, 1D heat equation (nonlinear)  
  -**Heat2D_pint_direct_fast.m**,  2D heat equation (linear)  
  -**Heat2D_pint_direct_nonlinear.m**, 2D heat equation (nonlinear)  
  
  -**Wave1D_pint_direct_fast.m**, 1D wave equation (linear)  
  -**fasteigB.m**, the core fast diagonalization algorithm  
  


  Licence Information:

This library is free software. It can be redistributed and/or modified under the terms of the MIT License.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Copyright (c) 2024 by Jun Liu
