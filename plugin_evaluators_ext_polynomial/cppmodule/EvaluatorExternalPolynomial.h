/*
Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
(HOOMD-blue) Open Source Software License Copyright 2009-2014 The Regents of
the University of Michigan All rights reserved.

HOOMD-blue may contain modifications ("Contributions") provided, and to which
copyright is held, by various Contributors who have granted The Regents of the
University of Michigan the right to modify and/or distribute such Contributions.

You may redistribute, use, and create derivate works of HOOMD-blue, in source
and binary forms, provided you abide by the following conditions:

* Redistributions of source code must retain the above copyright notice, this
list of conditions, and the following disclaimer both in the code and
prominently in any materials provided with the distribution.

* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions, and the following disclaimer in the documentation and/or
other materials provided with the distribution.

* All publications and presentations based on HOOMD-blue, including any reports
or published results obtained, in whole or in part, with HOOMD-blue, will
acknowledge its use according to the terms posted at the time of submission on:
http://codeblue.umich.edu/hoomd-blue/citations.html

* Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
http://codeblue.umich.edu/hoomd-blue/

* Apart from the above required attributions, neither the name of the copyright
holder nor the names of HOOMD-blue's contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id$
// $URL$
// Maintainer: joaander

#ifndef __EXTERNAL_EVALUATOR_POLYNOMIAL_H__
#define __EXTERNAL_EVALUATOR_POLYNOMIAL_H__

#ifndef NVCC
#include <string>
#endif

// always need to include hoomd_config first
#include <hoomd/hoomd_config.h>
#include <hoomd/HOOMDMath.h>

/*! \file EvaluatorExternalPolynomial.h
    \brief Defines the external evaluator class for Confinement potentials
    \details this also serves as the primary documetnation and
    base reference for the implementation of other external evaluators.
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
//! DEVICE is __host__ __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#else
#define DEVICE
#endif

//! Class for evaluating the confining external potential
/*! Note: This is identical to the pair potential EvaluatorPairLJ in hoomd. It is included here as an example.
    <b>General Overview</b>

    EvaluatorPairLJ is a low level computation class that computes the LJ pair potential V(r). As the standard
    MD potential, it also serves as a well documented example of how to write additional pair potentials. "Standard"
    pair potentials in hoomd are all handled via the template class PotentialPair. PotentialPair takes a potential
    evaluator as a template argument. In this way, all the complicated data mangament and other details of computing
    the pair force and potential on every single particle is only written once in the template class and the difference
    V(r) potentials that can be calculated are simply handled with various evaluator classes. Template instantiation is
    equivalent to inlining code, so there is no performance loss.

    In hoomd, a "standard" pair potential is defined as V(rsq, rcutsq, params, di, dj, qi, qj), where rsq is the squared
    distance between the two particles, rcutsq is the cuttoff radius at which the potential goes to 0, params is any
    number of per type-pair parameters, di, dj are the diameters of particles i and j, and qi, qj are the charges of
    particles i and j respectively.

    Diameter and charge are not always needed by a given pair evaluator, so it must provide the functions
    needsDiameter() and needsCharge() which return boolean values signifying if they need those quantities or not. A
    false return value notifies PotentialPair that it need not even load those valuse from memory, boosting performance.

    If needsDiameter() returns true, a setDiameter(Scalar di, Scalar dj) method will be called to set the two diameters.
    Similarly, if needsCharge() returns true, a setCharge(Scalar qi, Scalar qj) method will be called to set the two
    charges.

    All other arguments are common among all pair potentials and passed into the constructor. Coefficients are handled
    in a special way: the pair evaluator class (and PotentialPair) manage only a single parameter variable for each
    type pair. Pair potentials that need more than 1 parameter can specify that their param_type be a compound
    structure and reference that. For coalesced read performance on G200 GPUs, it is highly recommended that param_type
    is one of the following types: Scalar, Scalar2, Scalar4.

    The program flow will proceed like this: When a potential between a pair of particles is to be evaluated, a
    PairEvaluator is instantiated, passing the common parameters to the constructor and calling setDiameter() and/or
    setCharge() if need be. Then, the evalForceAndEnergy() method is called to evaluate the force and energy (more
    on that later). Thus, the evaluator must save all of the values it needs to compute the force and energy in member
    variables.

    evalForceAndEnergy() makes the necessary computations and sets the out parameters with the computed values.
    Specifically after the method complets, \a force_divr must be set to the value
    \f$ -\frac{1}{r}\frac{\partial V}{\partial r}\f$ and \a pair_eng must be set to the value \f$ V(r) \f$ if \a energy_shift is false or
    \f$ V(r) - V(r_{\mathrm{cut}}) \f$ if \a energy_shift is true.

    A pair potential evaluator class is also used on the GPU. So all of its members must be declared with the
    DEVICE keyword before them to mark them __device__ when compiling in nvcc and blank otherwise. If any other code
    needs to diverge between the host and device (i.e., to use a special math function like __powf on the device), it
    can similarly be put inside an ifdef NVCC block.

    <b>LJ specifics</b>

    EvaluatorPairLJ evaluates the function:
    \f[ V_{\mathrm{LJ}}(r) = 4 \varepsilon \left[ \left( \frac{\sigma}{r} \right)^{12} -
                                            \alpha \left( \frac{\sigma}{r} \right)^{6} \right] \f]
    broken up as follows for efficiency
    \f[ V_{\mathrm{LJ}}(r) = r^{-6} \cdot \left( 4 \varepsilon \sigma^{12} \cdot r^{-6} -
                                            4 \alpha \varepsilon \sigma^{6} \right) \f]
    . Similarly,
    \f[ -\frac{1}{r} \frac{\partial V_{\mathrm{LJ}}}{\partial r} = r^{-2} \cdot r^{-6} \cdot
            \left( 12 \cdot 4 \varepsilon \sigma^{12} \cdot r^{-6} - 6 \cdot 4 \alpha \varepsilon \sigma^{6} \right) \f]

    The LJ potential does not need diameter or charge. Two parameters are specified and stored in a Scalar2. \a lj1 is
    placed in \a params.x and \a lj2 is in \a params.y.

    These are related to the standard lj parameters sigma and epsilon by:
    - \a lj1 = 4.0 * epsilon * pow(sigma,12.0)
    - \a lj2 = alpha * 4.0 * epsilon * pow(sigma,6.0);

*/
class EvaluatorExternalPolynomial
    {
    public:
        //! Define the parameter type used by this pair potential evaluator
        typedef Scalar2 param_type;

        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance beteen the particles
            \param _rcutsq Sqauared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        DEVICE EvaluatorExternalPolynomial(Scalar3 X, const BoxDim& box, const param_type& params)
            : m_pos(X),
              m_box(box)
            {
            m_potentialcoeff= params.x;
            m_potentialpower = params.y;
            }

              //! Evaluate the force, energy and virial
        /*! \param F force vector
            \param energy value of the energy
            \param virial array of six scalars for the upper triangular virial tensor
        */
        DEVICE void evalForceEnergyAndVirial(Scalar3& F, Scalar& energy, Scalar* virial) {
            //Scalar3 a2 = make_scalar3(0,0,0);
            //Scalar3 a3 = make_scalar3(0,0,0);

            F.x = Scalar(0.0);
            F.y = Scalar(0.0);
            F.z = Scalar(0.0);
            energy = Scalar(0.0);

            // For this potential, since it uses scaled positions, the virial is always zero.
            for (unsigned int i = 0; i < 6; i++)
                virial[i] = Scalar(0.0);

            //Scalar V_box = m_box.getVolume();
         
            Scalar forcepower, forcecoeff, m_pos_squared, m_distance, fmag, costheta, sintheta;

            forcecoeff = m_potentialcoeff * m_potentialpower;
            forcepower = m_potentialpower - 1.0;
            m_pos_squared = dot(m_pos,m_pos);
            m_distance = fast::sqrt(m_pos_squared);
            costheta = m_pos.x / m_distance;
            sintheta = m_pos.y / m_distance;

            fmag = forcecoeff*fast::pow(m_distance, forcepower);
            //fmag = forcecoeff*fast::pow(m_pos.y,forcepower);

            F.x = Scalar(0);
            F.y = -1 * fmag * sintheta;
            //F.y = Scalar(-1)*Scalar(fmag);
            energy = m_potentialcoeff*fast::pow(m_pos.y,m_potentialpower);
        }

        #ifndef NVCC
        //! Get the name of this potential
        /*! \returns The potential name. Must be short and all lowercase, as this is the name energies will be logged as
            via analyze.log.
        */
        static std::string getName()
            {
            return std::string("polynomial");
            }
        #endif

    protected:
        Scalar3 m_pos;                //!< particle position
        BoxDim m_box;                 //!< box dimensions
        Scalar m_potentialcoeff;      //!< gamma
        Scalar m_potentialpower;      //!< n
    };


#endif // __EXTERNAL_EVALUATOR_POLYNOMIAL_H__
