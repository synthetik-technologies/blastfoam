/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "eigenSolver.H"
#include "scalar.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eigenSolver::eigenSolver
(
    const scalarSquareMatrix& A
)
:
    eigenvaluesRe_(A.n()),
    eigenvaluesIm_(A.n()),
    eigenvectors_(A.n(), A.n()),
    H_(),
    n_(A.n())
{
    if (isSymmetric(A))
    {
        eigenvectors_ = A;
        tridiagonaliseSymmMatrix();
        symmTridiagQL();
    }
    else
    {
        H_ = A;
        Hessenberg();
        realSchur();
    }
}

Foam::eigenSolver::eigenSolver(const scalarSquareMatrix& A, bool symmetric)
:
    eigenvaluesRe_(A.n()),
    eigenvaluesIm_(A.n()),
    eigenvectors_(A.n(), A.n()),
    H_(),
    n_(A.n())
{
    if (symmetric)
    {
        eigenvectors_ = A;
        tridiagonaliseSymmMatrix();
        symmTridiagQL();
    }
    else
    {
        H_ = A;
        Hessenberg();
        realSchur();
    }
}

bool Foam::eigenSolver::isSymmetric(const scalarSquareMatrix& A)
{
    bool symm = true;

    for (label j = 0; (j < n_) && symm; j++)
    {
        for (label i = 0; (i < n_) && symm; i++)
        {
            symm = (A[i][j] == A[j][i]);
        }
    }

    return symm;
}

void Foam::eigenSolver::tridiagonaliseSymmMatrix()
{
    for (label j = 0; j < n_; j++)
    {
        eigenvaluesRe_[j] = eigenvectors_[n_ - 1][j];
    }

    // Householder reduction to tridiagonal form
    for (label i = n_ - 1; i > 0; i--)
    {
        // Scale to avoid under/overflow.
        scalar scale = scalar(0);
        scalar h = scalar(0);

        for (label k = 0; k < i; k++)
        {
            scale = scale + mag(eigenvaluesRe_[k]);
        }

        if (scale == 0.0)
        {
            eigenvaluesIm_[i] = eigenvaluesRe_[i - 1];

            for (label j = 0; j < i; j++)
            {
                eigenvaluesRe_[j] = eigenvectors_[i - 1][j];
                eigenvectors_[i][j] = scalar(0);
                eigenvectors_[j][i] = scalar(0);
            }
        }
        else
        {
            // Generate Householder vector
            for (label k = 0; k < i; k++)
            {
                eigenvaluesRe_[k] /= scale;
                h += eigenvaluesRe_[k]*eigenvaluesRe_[k];
            }

            scalar f = eigenvaluesRe_[i - 1];
            scalar g = sqrt(h);

            if (f > 0)
            {
                g = -g;
            }

            eigenvaluesIm_[i] = scale*g;
            h -= f*g;
            eigenvaluesRe_[i - 1] = f - g;

            for (label j = 0; j < i; j++)
            {
                eigenvaluesIm_[j] = scalar(0);
            }

            // Apply similarity transformation to remaining columns
            for (label j = 0; j < i; j++)
            {
                f = eigenvaluesRe_[j];
                eigenvectors_[j][i] = f;
                g = eigenvaluesIm_[j] + eigenvectors_[j][j]*f;

                for (label k = j + 1; k <= i - 1; k++)
                {
                    g += eigenvectors_[k][j]*eigenvaluesRe_[k];
                    eigenvaluesIm_[k] += eigenvectors_[k][j]*f;
                }

                eigenvaluesIm_[j] = g;
            }

            f = scalar(0);

            for (label j = 0; j < i; j++)
            {
                eigenvaluesIm_[j] /= h;
                f += eigenvaluesIm_[j]*eigenvaluesRe_[j];
            }

            scalar hh = f/(2.0*h);

            for (label j = 0; j < i; j++)
            {
                eigenvaluesIm_[j] -= hh*eigenvaluesRe_[j];
            }

            for (label j = 0; j < i; j++)
            {
                f = eigenvaluesRe_[j];
                g = eigenvaluesIm_[j];

                for (label k = j; k <= i - 1; k++)
                {
                    eigenvectors_[k][j] -=
                            (f*eigenvaluesIm_[k] + g*eigenvaluesRe_[k]);
                }

                eigenvaluesRe_[j] = eigenvectors_[i - 1][j];
                eigenvectors_[i][j] = scalar(0);
            }
        }

        eigenvaluesRe_[i] = h;
    }

    // Accumulate transformations
    for (label i = 0; i < n_ - 1; i++)
    {
        eigenvectors_[n_ - 1][i] = eigenvectors_[i][i];
        eigenvectors_[i][i] = scalar(1);
        scalar h = eigenvaluesRe_[i + 1];

        if (h != 0.0)
        {
            for (label k = 0; k <= i; k++)
            {
                eigenvaluesRe_[k] = eigenvectors_[k][i + 1]/h;
            }

            for (label j = 0; j <= i; j++)
            {
                scalar g = scalar(0);

                for (label k = 0; k <= i; k++)
                {
                    g += eigenvectors_[k][i + 1]*eigenvectors_[k][j];
                }

                for (label k = 0; k <= i; k++)
                {
                    eigenvectors_[k][j] -= g*eigenvaluesRe_[k];
                }
            }
        }

        for (label k = 0; k <= i; k++)
        {
            eigenvectors_[k][i + 1] = scalar(0);
        }
    }

    for (label j = 0; j < n_; j++)
    {
        eigenvaluesRe_[j] = eigenvectors_[n_ - 1][j];
        eigenvectors_[n_ - 1][j] = scalar(0);
    }

    eigenvectors_[n_ - 1][n_ - 1] = scalar(1);
    eigenvaluesIm_[0] = scalar(0);
}

void Foam::eigenSolver::symmTridiagQL()
{
    for (label i = 1; i < n_; i++)
    {
        eigenvaluesIm_[i - 1] = eigenvaluesIm_[i];
    }

    eigenvaluesIm_[n_-1] = scalar(0);

    scalar f = scalar(0);
    scalar tst1 = scalar(0);
    scalar eps = pow(2.0, -52.0);

    for (label l = 0; l < n_; l++)
    {
        // Find small subdiagonal element
        tst1 = max(tst1, mag(eigenvaluesRe_[l]) + mag(eigenvaluesIm_[l]));
        label m = l;

        // Original while-loop from Java code
        while (m < n_)
        {
            if (mag(eigenvaluesIm_[m]) <= eps*tst1)
            {
                break;
            }

            m++;
        }

        // If m == l, eigenvaluesRe_[l] is an eigenvalue, otherwise, iterate
        if (m > l)
        {
            label iter = 0;

            do
            {
                iter += 1;

                // Compute implicit shift
                scalar g = eigenvaluesRe_[l];
                scalar p = (eigenvaluesRe_[l + 1] - g)/(2.0*eigenvaluesIm_[l]);
                scalar r = hypot(p, 1.0);

                if (p < 0)
                {
                    r = -r;
                }

                eigenvaluesRe_[l] = eigenvaluesIm_[l]/(p + r);
                eigenvaluesRe_[l + 1] = eigenvaluesIm_[l]*(p + r);
                scalar dl1 = eigenvaluesRe_[l + 1];
                scalar h = g - eigenvaluesRe_[l];

                for (label i = l + 2; i < n_; i++)
                {
                    eigenvaluesRe_[i] -= h;
                }

                f += h;

                // Implicit QL transformation.
                p = eigenvaluesRe_[m];
                scalar c = scalar(1);
                scalar c2 = c;
                scalar c3 = c;
                scalar el1 = eigenvaluesIm_[l + 1];
                scalar s = scalar(0);
                scalar s2 = scalar(0);

                for (label i = m - 1; i >= l; i--)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c*eigenvaluesIm_[i];
                    h = c*p;
                    r = hypot(p, eigenvaluesIm_[i]);
                    eigenvaluesIm_[i + 1] = s*r;
                    s = eigenvaluesIm_[i]/r;
                    c = p/r;
                    p = c*eigenvaluesRe_[i] - s*g;
                    eigenvaluesRe_[i + 1] = h + s*(c*g + s*eigenvaluesRe_[i]);

                    // Accumulate transformation
                    for (label k = 0; k < n_; k++)
                    {
                        h = eigenvectors_[k][i + 1];
                        eigenvectors_[k][i + 1] = s*eigenvectors_[k][i] + c*h;
                        eigenvectors_[k][i] = c*eigenvectors_[k][i] - s*h;
                    }
                }

                p = -s*s2*c3*el1*eigenvaluesIm_[l]/dl1;
                eigenvaluesIm_[l] = s*p;
                eigenvaluesRe_[l] = c*p;
            }
            while (mag(eigenvaluesIm_[l]) > eps*tst1); // Convergence check
        }

        eigenvaluesRe_[l] = eigenvaluesRe_[l] + f;
        eigenvaluesIm_[l] = scalar(0);
    }

    // Sorting eigenvalues and corresponding vectors
    for (label i = 0; i < n_ - 1; i++)
    {
        label k = i;
        scalar p = eigenvaluesRe_[i];

        for (label j = i + 1; j < n_; j++)
        {
            if (eigenvaluesRe_[j] < p)
            {
                k = j;
                p = eigenvaluesRe_[j];
            }
        }

        if (k != i)
        {
            eigenvaluesRe_[k] = eigenvaluesRe_[i];
            eigenvaluesRe_[i] = p;

            for (label j = 0; j < n_; j++)
            {
                p = eigenvectors_[j][i];
                eigenvectors_[j][i] = eigenvectors_[j][k];
                eigenvectors_[j][k] = p;
            }
        }
    }
}

void Foam::eigenSolver::Hessenberg()
{
    scalarDiagonalMatrix orth_(n_);

    label low = 0;
    label high = n_ - 1;

    for (label m = low + 1; m <= high - 1; m++)
    {
        // Scale column
        scalar scale = scalar(0);

        for (label i = m; i <= high; i++)
        {
            scale = scale + mag(H_[i][m - 1]);
        }

        if (scale != 0.0)
        {
            // Compute Householder transformation
            scalar h = scalar(0);

            for (label i = high; i >= m; i--)
            {
                orth_[i] = H_[i][m - 1]/scale;
                h += orth_[i]*orth_[i];
            }

            scalar g = sqrt(h);

            if (orth_[m] > 0)
            {
                g = -g;
            }

            h -= orth_[m]*g;
            orth_[m] -= g;

            // Apply Householder similarity transformation
            // H = (I - u*u'/h)*H*(I - u*u')/h)
            for (label j = m; j < n_; j++)
            {
                scalar f = scalar(0);

                for (label i = high; i >= m; i--)
                {
                    f += orth_[i]*H_[i][j];
                }

                f /= h;

                for (label i = m; i <= high; i++)
                {
                    H_[i][j] -= f*orth_[i];
                }
            }

            for (label i = 0; i <= high; i++)
            {
                scalar f = scalar(0);

                for (label j = high; j >= m; j--)
                {
                    f += orth_[j]*H_[i][j];
                }

                f /= h;

                for (label j = m; j <= high; j++)
                {
                    H_[i][j] -= f*orth_[j];
                }
            }

            orth_[m] = scale*orth_[m];
            H_[m][m-1] = scale*g;
        }
    }

    // Accumulate transformations
    for (label i = 0; i < n_; i++)
    {
        for (label j = 0; j < n_; j++)
        {
            eigenvectors_[i][j] = (i == j ? 1.0 : 0.0);
        }
    }

    for (label m = high - 1; m >= low + 1; m--)
    {
        if (H_[m][m - 1] != 0.0)
        {
            for (label i = m + 1; i <= high; i++)
            {
                orth_[i] = H_[i][m - 1];
            }

            for (label j = m; j <= high; j++)
            {
                scalar g = scalar(0);

                for (label i = m; i <= high; i++)
                {
                    g += orth_[i]*eigenvectors_[i][j];
                }

                // Double division avoids possible underflow
                g = (g/orth_[m])/H_[m][m - 1];

                for (label i = m; i <= high; i++)
                {
                    eigenvectors_[i][j] += g*orth_[i];
                }
            }
        }
    }
}

void Foam::eigenSolver::realSchur()
{
    // Initialize
    label nn = n_;
    label n = nn - 1;
    label low = 0;
    label high = nn - 1;
    scalar eps = pow(2.0, -52.0);
    scalar exshift = scalar(0);
    scalar p = 0, q = 0, r = 0, s = 0, z = 0, t, w, x, y;

    // Store roots isolated by balance and compute matrix norm
    scalar norm = scalar(0);

    for (label i = 0; i < nn; i++)
    {
        if ((i < low) || (i > high))
        {
            eigenvaluesRe_[i] = H_[i][i];
            eigenvaluesIm_[i] = scalar(0);
        }

        for (label j = max(i - 1, 0); j < nn; j++)
        {
            norm += mag(H_[i][j]);
        }
    }

    label iter = 0;

    while (n >= low)
    {
        // Look for single small sub-diagonal element
        label l = n;

        while (l > low)
        {
            s = mag(H_[l - 1][l - 1]) + mag(H_[l][l]);

            if (s == 0.0)
            {
                s = norm;
            }

            if (mag(H_[l][l - 1]) < eps*s)
            {
                break;
            }

            l--;
        }

        // Check for convergence
        if (l == n)         // One root found
        {
            H_[n][n] = H_[n][n] + exshift;
            eigenvaluesRe_[n] = H_[n][n];
            eigenvaluesIm_[n] = scalar(0);
            n--;
            iter = 0;
        }
        else if (l == n - 1)  // Two roots found
        {
            w = H_[n][n - 1]*H_[n - 1][n];
            p = (H_[n - 1][n - 1] - H_[n][n])/2.0;
            q = p*p + w;
            z = sqrt(mag(q));
            H_[n][n] += exshift;
            H_[n - 1][n - 1] += exshift;
            x = H_[n][n];

            // Scalar pair
            if (q >= 0)
            {
                if (p >= 0)
                {
                    z = p + z;
                }
                else
                {
                    z = p - z;
                }

                eigenvaluesRe_[n - 1] = x + z;
                eigenvaluesRe_[n] = eigenvaluesRe_[n - 1];

                if (z != 0.0)
                {
                    eigenvaluesRe_[n] = x - w/z;
                }

                eigenvaluesIm_[n - 1] = scalar(0);
                eigenvaluesIm_[n] = scalar(0);
                x = H_[n][n - 1];
                s = mag(x) + mag(z);
                p = x/s;
                q = z/s;
                r = sqrt(p*p + q*q);
                p /= r;
                q /= r;

                // Row modification
                for (label j = n - 1; j < nn; j++)
                {
                    z = H_[n - 1][j];
                    H_[n - 1][j] = q*z + p*H_[n][j];
                    H_[n][j] = q*H_[n][j] - p*z;
                }

                // Column modification
                for (label i = 0; i <= n; i++)
                {
                    z = H_[i][n - 1];
                    H_[i][n - 1] = q*z + p*H_[i][n];
                    H_[i][n] = q*H_[i][n] - p*z;
                }

                // Accumulate transformations
                for (label i = low; i <= high; i++)
                {
                    z = eigenvectors_[i][n - 1];
                    eigenvectors_[i][n - 1] = q*z + p*eigenvectors_[i][n];
                    eigenvectors_[i][n] = q*eigenvectors_[i][n] - p*z;
                }
            }
            else         // Complex pair of roots
            {
                eigenvaluesRe_[n - 1] = x + p;
                eigenvaluesRe_[n] = x + p;
                eigenvaluesIm_[n - 1] = z;
                eigenvaluesIm_[n] = -z;
            }

            n -= 2;
            iter = 0;
        }
        else            // No convergence yet
        {
            // Form shift
            x = H_[n][n];
            y = scalar(0);
            w = scalar(0);

            if (l < n)
            {
                y = H_[n - 1][n - 1];
                w = H_[n][n - 1]*H_[n - 1][n];
            }

            // Wilkinson's original ad-hoc shift
            if (iter == 10)
            {
                exshift += x;
                for (label i = low; i <= n; i++)
                {
                    H_[i][i] -= x;
                }

                s = mag(H_[n][n - 1]) + mag(H_[n - 1][n - 2]);
                x = 0.75*s;
                y = 0.75*s;
                w = -0.4375*s*s;
            }

            // New ad hoc shift
            if (iter == 30)
            {
                s = (y - x)/2.0;
                s = s*s + w;

                if (s > 0)
                {
                    s = sqrt(s);

                    if (y < x)
                    {
                        s = -s;
                    }

                    s = x - w/((y - x)/2.0 + s);

                    for (label i = low; i <= n; i++)
                    {
                        H_[i][i] -= s;
                    }

                    exshift += s;
                    x = 0.964;
                    y = 0.964;
                    w = 0.964;
                }
            }

            iter += 1;

            // Look for two consecutive small sub-diagonal elements
            label m = n - 2;

            while (m >= l)
            {
                z = H_[m][m];
                r = x - z;
                s = y - z;
                p = (r*s - w)/H_[m + 1][m] + H_[m][m + 1];
                q = H_[m + 1][m + 1] - z - r - s;
                r = H_[m + 2][m + 1];
                s = mag(p) + mag(q) + mag(r);
                p /= s;
                q /= s;
                r /= s;

                if (m == l)
                {
                    break;
                }

                if
                (
                    mag(H_[m][m - 1])*(mag(q) + mag(r))
                            < eps*(mag(p)*(mag(H_[m - 1][m - 1])
                                    + mag(z) + mag(H_[m + 1][m + 1])))
                )
                {
                    break;
                }

                m--;
            }

            for (label i = m + 2; i <= n; i++)
            {
                H_[i][i - 2] = scalar(0);

                if (i > m + 2)
                {
                    H_[i][i - 3] = scalar(0);
                }
            }

            // Double QR step involving rows l:n and columns m:n
            for (label k = m; k <= n - 1; k++)
            {
                label notlast = (k != n - 1);

                if (k != m)
                {
                    p = H_[k][k - 1];
                    q = H_[k + 1][k - 1];
                    r = (notlast ? H_[k + 2][k - 1] : 0.0);
                    x = mag(p) + mag(q) + mag(r);

                    if (x != 0.0)
                    {
                        p /= x;
                        q /= x;
                        r /= x;
                    }
                }

                if (x == 0.0)
                {
                    break;
                }

                s = sqrt(p*p + q*q + r*r);

                if (p < 0)
                {
                    s = -s;
                }

                if (s != 0)
                {
                    if (k != m)
                    {
                        H_[k][k - 1] = -s*x;
                    }
                    else if (l != m)
                    {
                        H_[k][k - 1] = -H_[k][k - 1];
                    }

                    p = p + s;
                    x = p/s;
                    y = q/s;
                    z = r/s;
                    q /= p;
                    r /= p;

                    // Row modification
                    for (label j = k; j < nn; j++)
                    {
                        p = H_[k][j] + q*H_[k + 1][j];

                        if (notlast)
                        {
                            p += r*H_[k + 2][j];
                            H_[k + 2][j] -= p*z;
                        }

                        H_[k][j] -= p*x;
                        H_[k + 1][j] -= p*y;
                    }

                    // Column modification
                    for (label i = 0; i <= min(n, k + 3); i++)
                    {
                        p = x*H_[i][k] + y*H_[i][k + 1];

                        if (notlast)
                        {
                            p += z*H_[i][k + 2];
                            H_[i][k + 2] -= p*r;
                        }

                        H_[i][k] -= p;
                        H_[i][k + 1] -= p*q;
                    }

                    // Accumulate transformations
                    for (label i = low; i <= high; i++)
                    {
                        p = x*eigenvectors_[i][k] + y*eigenvectors_[i][k + 1];

                        if (notlast)
                        {
                            p += z*eigenvectors_[i][k + 2];
                            eigenvectors_[i][k + 2] -= p*r;
                        }

                        eigenvectors_[i][k] -= p;
                        eigenvectors_[i][k + 1] -= p*q;
                    }
                }
            }
        }
    }

    // Backsubstitute to find vectors of upper triangular form
    if (norm == 0.0)
    {
        return;
    }

    for (n = nn-1; n >= 0; n--)
    {
        p = eigenvaluesRe_[n];
        q = eigenvaluesIm_[n];

        // scalar vector
        if (q == 0)
        {
            label l = n;
            H_[n][n] = scalar(1);

            for (label i = n-1; i >= 0; i--)
            {
                w = H_[i][i] - p;
                r = scalar(0);

                for (label j = l; j <= n; j++)
                {
                    r += H_[i][j]*H_[j][n];
                }

                if (eigenvaluesIm_[i] < 0.0)
                {
                    z = w;
                    s = r;
                }
                else
                {
                    l = i;

                    if (eigenvaluesIm_[i] == 0.0)
                    {
                        if (w != 0.0)
                        {
                            H_[i][n] = -r/w;
                        }
                        else
                        {
                            H_[i][n] = -r/(eps*norm);
                        }
                    }
                    else // Solve real equations
                    {
                        x = H_[i][i + 1];
                        y = H_[i + 1][i];
                        q = (eigenvaluesRe_[i] - p)*(eigenvaluesRe_[i] - p)
                                + eigenvaluesIm_[i]*eigenvaluesIm_[i];

                        t = (x*s - z*r)/q;
                        H_[i][n] = t;

                        if (mag(x) > mag(z))
                        {
                            H_[i + 1][n] = (-r - w*t)/x;
                        }
                        else
                        {
                            H_[i + 1][n] = (-s - y*t)/z;
                        }
                    }

                    // Overflow control
                    t = mag(H_[i][n]);

                    if ((eps*t)*t > 1)
                    {
                        for (label j = i; j <= n; j++)
                        {
                            H_[j][n] /= t;
                        }
                    }
                }
            }
        }
        else if (q < 0) // Complex vector
        {
            label l = n - 1;

            // Last vector component imaginary so matrix is triangular
            if (mag(H_[n][n - 1]) > mag(H_[n - 1][n]))
            {
                H_[n - 1][n - 1] = q/H_[n][n - 1];
                H_[n - 1][n] = -(H_[n][n] - p)/H_[n][n - 1];
            }
            else
            {
                complex cDiv = complex(0.0, -H_[n - 1][n])
                        /complex(H_[n - 1][n-1]-p, q);

                H_[n - 1][n - 1] = cDiv.Re();
                H_[n - 1][n] = cDiv.Im();
            }

            H_[n][n - 1] = scalar(0);
            H_[n][n] = scalar(1);

            for (label i = n - 2; i >= 0; i--)
            {
                scalar ra, sa, vr, vi;
                ra = scalar(0);
                sa = scalar(0);

                for (label j = l; j <= n; j++)
                {
                    ra += H_[i][j]*H_[j][n-1];
                    sa += H_[i][j]*H_[j][n];
                }

                w = H_[i][i] - p;

                if (eigenvaluesIm_[i] < 0.0)
                {
                    z = w;
                    r = ra;
                    s = sa;
                }
                else
                {
                    l = i;

                    if (eigenvaluesIm_[i] == 0)
                    {
                        complex cDiv = complex(-ra, -sa)/complex(w, q);
                        H_[i][n - 1] = cDiv.Re();
                        H_[i][n] = cDiv.Im();
                    }
                    else
                    {
                        // Solve complex equations
                        x = H_[i][i + 1];
                        y = H_[i + 1][i];
                        vr = (eigenvaluesRe_[i] - p)*(eigenvaluesRe_[i] - p)
                                + eigenvaluesIm_[i]*eigenvaluesIm_[i] - q*q;

                        vi = 2.0*(eigenvaluesRe_[i] - p)*q;

                        if ((vr == 0.0) && (vi == 0.0))
                        {
                            vr = eps*norm*(mag(w) + mag(q) + mag(x) + mag(y)
                                    + mag(z));
                        }

                        complex cDiv =
                                complex(x*r - z*ra + q*sa, x*s - z*sa - q*ra)
                                    /complex(vr, vi);

                        H_[i][n - 1] = cDiv.Re();
                        H_[i][n] = cDiv.Im();

                        if (mag(x) > (mag(z) + mag(q)))
                        {
                            H_[i + 1][n - 1] = (-ra - w*H_[i][n - 1]
                                    + q*H_[i][n])/x;

                            H_[i + 1][n] = (-sa - w*H_[i][n]
                                    - q*H_[i][n - 1])/x;
                        }
                        else
                        {
                            complex cDiv = complex(-r - y*H_[i][n - 1], -s
                                    - y*H_[i][n])/complex(z, q);

                            H_[i + 1][n - 1] = cDiv.Re();
                            H_[i + 1][n] = cDiv.Im();
                        }
                    }

                    // Overflow control
                    t = max(mag(H_[i][n - 1]), mag(H_[i][n]));

                    if ((eps*t)*t > 1)
                    {
                        for (label j = i; j <= n; j++)
                        {
                            H_[j][n - 1] /= t;
                            H_[j][n] /= t;
                        }
                    }
                }
            }
        }
    }

    // Vectors of isolated roots
    for (label i = 0; i < nn; i++)
    {
        if (i < low || i > high)
        {
            for (label j = i; j < nn; j++)
            {
                eigenvectors_[i][j] = H_[i][j];
            }
        }
    }

    // Back transformation to get eigenvectors of original matrix
    for (label j = nn - 1; j >= low; j--)
    {
        for (label i = low; i <= high; i++)
        {
            z = scalar(0);

            for (label k = low; k <= min(j, high); k++)
            {
                z = z + eigenvectors_[i][k]*H_[k][j];
            }

            eigenvectors_[i][j] = z;
        }
    }
}

// ************************************************************************* //
