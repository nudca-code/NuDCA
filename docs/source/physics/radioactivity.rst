=================
Radioactive Decay
=================

Radioactive decay chains describe the sequential transformation of unstable nuclei into more stable ones through a series of radioactive decays. Each nuclide in the chain decays into the next, until a stable nuclide is reached. The mathematical description of these chains is essential for understanding the time evolution of the number of nuclei of each type.

**The decay chain equations**

A decay process with branching decay chains can be described as:


.. math::
    \begin{cases}
        \frac{\mathrm{d}N_{1}(t)}{\mathrm{d}t} &= -\lambda_{1}N_{1}(t) \\
        \frac{\mathrm{d}N_{2}(t)}{\mathrm{d}t} &= b_{12}\lambda_{1}N_{1}(t) - \lambda_{2}N_{2}(t) \\
        \frac{\mathrm{d}N_{3}(t)}{\mathrm{d}t} &= b_{13}\lambda_{1}N_{1}(t) + b_{23}\lambda_{2}N_{2}(t) - \lambda_{3}N_{3}(t) \\
        &\vdots \\
        \frac{\mathrm{d}N_{n}(t)}{\mathrm{d}t} &= b_{1n}\lambda_{1}N_{1}(t) + b_{2n}\lambda_{2}N_{2}(t) + \cdots + b_{n-1,n}\lambda_{n-1}N_{n-1}(t) - \lambda_{n}N_{n}(t)
    \end{cases}

Here:

- :math:`N_i(t)` is the number of nuclei of type :math:`i` at time :math:`t`.
- :math:`\lambda_i` is the decay constant for nuclide :math:`i`.
- :math:`b_{ji}` is the branching ratio for the decay from nuclide :math:`j` to nuclide :math:`i`.


**Matrix Representation**

The above system can be written in matrix form as:

.. math::
    \frac{\mathrm{d}}{\mathrm{d}t}
        \begin{bmatrix}
            N_{1}(t) \\
            N_{2}(t) \\
            N_{3}(t) \\
            \vdots   \\
            N_{i}(t) \\
            \vdots   \\
            N_{n}(t)
        \end{bmatrix} =
        \begin{bmatrix}
            -\lambda_{1}      &                   &                   &          &                    &        &               \\
            b_{12}\lambda_{1} & -\lambda_{2}      &                   &          &                    &        &               \\
            b_{13}\lambda_{1} & b_{23}\lambda_{2} & -\lambda_{3}      &          &                    &        &               \\
            \vdots            & \vdots            & \vdots            & \ddots   &                    &        &               \\
            b_{1i}\lambda_{1} & b_{2i}\lambda_{2} & b_{3i}\lambda_{3} & \cdots   & -\lambda_{i}       &        &               \\
            \vdots            & \vdots            & \vdots            & \vdots   & \vdots             & \ddots &               \\
            b_{1n}\lambda_{1} & b_{2n}\lambda_{2} & b_{3n}\lambda_{3} & \cdots   &  b_{in}\lambda_{i} & \cdots & -\lambda_{n}
        \end{bmatrix}
        \times
        \begin{bmatrix}
            N_{1}(t) \\
            N_{2}(t) \\
            N_{3}(t) \\
            \vdots   \\
            N_{i}(t) \\
            \vdots   \\
            N_{n}(t)
        \end{bmatrix}


where :math:`\boldsymbol{A}` is the transition matrix:

.. math::
    \boldsymbol{A} = 
    \begin{bmatrix}
        -\lambda_{1}      &                   &                   &          &                    &        &                \\
        b_{12}\lambda_{1} & -\lambda_{2}      &                   &          &                    &        &                \\
        b_{13}\lambda_{1} & b_{23}\lambda_{2} & -\lambda_{3}      &          &                    &        &                \\
        \vdots            & \vdots            & \vdots            & \ddots   &                    &        &                \\
        b_{1i}\lambda_{1} & b_{2i}\lambda_{2} & b_{3i}\lambda_{3} & \cdots   & -\lambda_{i}       &        &                \\
        \vdots            & \vdots            & \vdots            & \vdots   & \vdots             & \ddots &                \\
        b_{1n}\lambda_{1} & b_{2n}\lambda_{2} & b_{3n}\lambda_{3} & \cdots   &  b_{in}\lambda_{i} & \cdots & -\lambda_{n}
    \end{bmatrix}


The elements of the matrix :math:`A_{ij}` are defined as:

.. math::
    A_{ij} =
    \begin{cases} 
        0                 & \text{if } i < j, \\
        -\lambda_{i}      & \text{if } i = j, \\
        b_{ji}\lambda_{j} & \text{if } i > j.
    \end{cases}

This lower-triangular structure reflects the fact that each nuclide can only be produced by the decay of its precursors.

**General Solution**

The system can be compactly written as:

.. math::
    \frac{\mathrm{d}\boldsymbol{N}}{\mathrm{d}t} = \boldsymbol{A} \boldsymbol{N}



The general solution involves diagonalizing the matrix :math:`\boldsymbol{A}`, which requires finding its eigenvalues and eigenvectors.

**Eigenvalues and Eigenvectors**

The eigenvalues :math:`\alpha_j` of :math:`\boldsymbol{A}` are given by:

.. math::
    \alpha_{j} = -\lambda_{j} = A_{jj} \quad(j = 1,2,3,\dots,n)

The corresponding eigenvectors :math:`\boldsymbol{p}_j` have the form:

.. math::
    \boldsymbol{p}_{j} = 
    \begin{bmatrix}
        P_{1j}    \\
        P_{2j}    \\
        P_{3j}    \\
        \vdots   \\
        P_{ij}    \\
        \vdots   \\
        P_{nj}
    \end{bmatrix}

where the components :math:`P_{ij}` satisfy the following recurrence relations:

.. math::
    \begin{cases}
        (A_{11} - A_{jj})P_{1j} &= 0                           \\
        A_{21}P_{1j} + (A_{22} - A_{jj})P_{2j} &= 0        \\
        A_{31}P_{1j} + A_{32}P_{2j} + (A_{33} - A_{jj})P_{3j} &= 0         \\
        \vdots                                           \\
        A_{n1}P_{1j} + A_{n2}P_{2j} + A_{n3}P_{3j} + \cdots + (A_{nn} - A_{jj})P_{nj} &= 0
    \end{cases}


Or, in matrix form:

.. math::
    \begin{bmatrix}
        A_{11} - A_{jj} &                 &                 &              &                 &                    &                                                   \\
        A_{21}          & A_{22} - A_{jj} &                 &              &                        &                    &                     \\
        A_{31}          & A_{32}          & A_{33} - A_{jj} &              &                        &               &                     \\
        \vdots            & \vdots            & \vdots            & \ddots       &                  &                    &        &              \\
        A_{i1}          & A_{i2}          & A_{i3}          & \cdots   & A_{ii} - A_{jj}       &        &               \\
        \vdots            & \vdots            & \vdots            & \vdots       & \vdots            & \ddots             &                \\
        A_{n1} & A_{n2} & A_{n3} & \cdots       &  A_{ni} & \cdots & A_{nn} - A_{jj}
    \end{bmatrix}
    \begin{bmatrix}
        P_{1j}    \\
        P_{2j}    \\
        P_{3j}    \\
        \vdots   \\
        P_{ij}    \\
        \vdots   \\
        P_{nj}
    \end{bmatrix} = 
    \begin{bmatrix}
        0   \\
        0   \\
        0   \\
        \vdots \\
        0    \\
        \vdots  \\
        0
    \end{bmatrix}

The explicit solution for :math:`P_{ij}` is:

.. math::
    P_{ij} =
    \begin{cases} 
        0 & \text{if } i < j, \\
        1 & \text{if } i = j, \\
        \frac{1}{A_{jj} - A_{ii}} \sum_{k=j}^{i-1} A_{ik}P_{kj} & \text{if } i > j.
    \end{cases}


Similarly, the inverse matrix elements :math:`P^{-1}_{ij}` are given by:

.. math::
    P^{-1}_{ij} =
    \begin{cases} 
        0 & \text{if } i < j, \\
        1 & \text{if } i = j, \\
        -\sum_{k=j}^{i-1} P_{ik}P^{-1}_{kj} & \text{if } i > j.
    \end{cases}


This formalism provides a general framework for analyzing arbitrary radioactive decay chains, including those with branching decays. The matrix approach is particularly powerful for numerical computation and for understanding the structure of the solutions.
