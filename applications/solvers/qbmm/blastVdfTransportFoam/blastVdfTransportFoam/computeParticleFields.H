labelList orderZero(momentOrders[0].size(), 0);
volScalarField m0(max(moments(orderZero), 1e-10));

forAll(velocityIndexes, cmpt)
{
    labelList orderOne(orderZero);
    orderOne[velocityIndexes[cmpt]] = 1;

    volScalarField meanU(moments(orderOne)/m0);
    Up.replace(cmpt, meanU);

    if (computeVariance)
    {
        forAll(velocityIndexes, cmpt2)
        {
            labelList orderOne2(orderZero);
            labelList orderTwo(orderZero);
            orderOne2[velocityIndexes[cmpt2]] = 1;
            orderTwo[velocityIndexes[cmpt]] = 1;
            orderTwo[velocityIndexes[cmpt2]] += 1;

            volScalarField meanU2(moments(orderOne2)/m0);

            volScalarField coVar
            (
                moments(orderTwo)/m0 - meanU*meanU2
            );
            Sigmap().replace(cmpt + cmpt2*3, coVar);
        }
    }
}

if (computeVariance)
{
    Thetap() = tr(Sigmap())/3.0;
}
