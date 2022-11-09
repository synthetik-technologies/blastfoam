#include "../../../MacDonald.H"
scalar S0Eqn::h(const scalar x) const
{
    static const scalar a1 = 0.674202;
    static const scalar a2 = 21.7112;
    static const scalar a3 = 14.492;
    static const scalar a4 = 1.4305;
    scalar xBar = x/100.0 - 2.0/3.0;
    return
        (x < (200.0/3.0))
      ? (
            cbrt(4.0/g_)*(4.0/3.0 - x/100.0)
          - 9.0*x/1000.0*xBar
        )
      : cbrt(4.0/g_)
       *(
            a1*(pow4(xBar) + pow3(xBar))
          - a2*sqr(xBar)
          + a3*xBar
          + a4
        );
}
scalar S0Eqn::dhdx(const scalar x) const
{
    static const scalar a1 = 0.674202;
    static const scalar a2 = 21.7112;
    static const scalar a3 = 14.492;
    static const scalar a4 = 1.4305;
    scalar xBar = x/100.0 - 2.0/3.0;
    return
        (x < (200.0/3.0))
      ? (
          - cbrt(4.0/g_)/100.0
          - 9.0/500.0*(x/100.0 - 1.0/3.0)
        )
      : cbrt(4.0/g_)
       *(
            a1*(4.0*pow3(xBar) + 3.0*sqr(xBar))
          - a2*2.0*xBar
          + a3
        )/100.0;
}