#include "IntegratorsFwd.H"
#include "createEquations.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "OFstream.H"

namespace Foam
{

class S0Eqn
:
    public ScalarEquation
{
    const fileName caseDir_;
    const dictionary dict_;
    const scalar g_;
    const scalar L_;
    const scalar q_;
    scalar n_;
    scalar f_;
    scalar R_;
    scalar muh_;
    scalar muv_;
    scalar kl_;
    scalar kt_;

public:

    S0Eqn()
    :
        caseDir_(getEnv("FOAM_CASE")),
        dict_((IFstream(caseDir_/"constant/shallowWaterProperties"))()),
        g_(mag(dict_.lookup<vector>("g"))),
        L_
        (
            dictionary
            (
                (IFstream(caseDir_/"system/blockMeshDict"))()
            ).lookup<List<vector>>("vertices")[1].x()
        ),
        q_
        (
            vectorField
            (
                "inletValue",
                dictionary
                (
                    (IFstream(caseDir_/"0/U"))()
                ).subDict
                (
                    "boundaryField"
                ).subDict("inlet"),
                1
            )[0].x()
        ),
        n_(0.0),
        f_(0.0),

        R_(0.0),

        muh_(0.0),
        muv_(0.0),
        kl_(0.0),
        kt_(0.0)


    {
        const bool friction = dict_.lookupOrDefault("friction", false);
        const bool rain = dict_.lookupOrDefault("rain", false);
        const bool viscous = dict_.lookupOrDefault("viscous", false);
        if (friction || rain || viscous)
        {
            if (friction)
            {
                word fType(dict_.lookup("frictionType"));
                if (fType == "Manning")
                {
                    n_ = dict_.lookup<scalar>("n");
                }
                else
                {
                    f_ = dict_.lookup<scalar>("f");
                }
            }
            if (rain)
            {
                R_ = dict_.lookup<scalar>("R");
            }
            if (viscous)
            {
                const bool tubulence = dict_.lookup<bool>("turbulence");
                if (tubulence)
                {
                    muh_ = dict_.lookup<scalar>("muh");
                    muv_ = dict_.lookup<scalar>("muv");
                    kl_ = dict_.lookup<scalar>("kl");
                    kt_ = dict_.lookup<scalar>("kt");
                }
                else
                {
                    muh_ = dict_.lookup<scalar>("mu")/4.0;
                }
            }
        }
    }

    label nDerivatives() const
    {
        return 0;
    }
    scalar L() const
    {
        return L_;
    }
    scalar q(const scalar x) const
    {
        return q_ + x*R_;
    }

    scalar h(const scalar x) const;
    scalar dhdx(const scalar x) const;
    scalar d2hdx2(const scalar x) const;

    scalar fx(const scalar x, const label li) const
    {
        scalar S0 = (sqr(q(x))/(g_*pow3(h(x))) - 1.0)*dhdx(x);
        if (n_ > 0)
        {
            S0 -= sqr(n_)*sqr(q(x))/pow(h(x), 10.0/3.0);
        }
        if (f_ > 0)
        {
            S0 -= f_*sqr(q(x))/(8.0*g_*pow3(h(x)));
        }
        if (R_ > 0)
        {
            S0 -= 2.0*q(x)*R_/(g_*sqr(h(x)));
        }
        if (muh_ > 0)
        {
            S0 -= 4.0*muh_/(g_*sqr(h(x)))*(q(x)*d2hdx2(x) + q(x)/h(x)*sqr(dhdx(x)));
        }

        if (muv_ > 0)
        {
            if (kl_ > 0)
            {
                S0 -= kl_*q(x)/(g_*sqr(h(x))*(1.0 + kl_*h(x)/(3.0*muv_)));
                if (kt_ > 0)
                {
                    S0 -= kt_*mag(q(x))*q(x)/(g_*sqr(h(x))*sqr(1.0 + kl_*h(x)/(3.0*muv_)));
                }
            }
        }
        return S0;
    }

    scalar h0(const scalar x, const bool adjust = true) const
    {
        static autoPtr<scalarIntegrator> integrator;
        if (!integrator.valid())
        {
            integrator =
                scalarIntegrator::New("Simpson38", *this, dictionary());
            OFstream hos(getEnv("FOAM_CASE")/"h");
            hos << "h_L uniform " << h(L_) << ';' << endl;
            hos << "h_0 uniform " << h(0) << ';' << endl;
        }
        return
            (adjust ? integrator->integrate(0, L_) : 0.0)
          - integrator->integrate(x, L_);
    }
};

#ifndef d2dx2
scalar S0Eqn::d2hdx2(const scalar x) const
{
    return 0;
}
#endif
}

using namespace Foam;
// ************************************************************************* //
