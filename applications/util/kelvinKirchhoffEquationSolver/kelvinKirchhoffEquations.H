class kelvinKirchhoffEquations
:
    public ODESystem
{
    //- Handle to floaterMotion whose calcAcceleration function is used and from
    //  whose input file data is read.
    floaterMotion& body_;
    const labelList DoFs_;

public:

    kelvinKirchhoffEquations(floaterMotion& body, const labelList& DoFs)
    :
    body_(body),
    DoFs_(DoFs)
    {}

    label nEqns() const
    {
        return 18;
    }

    void derivatives
    (
        const scalar t,
        const scalarField& Y,
        scalarField& dYdt
    ) const
    {

        // State: x0, Q, v0, omega
        vector x0(Y[0], Y[1], Y[2]);
        tensor Q
        (
            Y[3], Y[4], Y[5],
            Y[6], Y[7], Y[8],
            Y[9], Y[10], Y[11]
        );
        vector v0(Y[12], Y[13], Y[14]);
        vector omega(Y[15], Y[16], Y[17]);

        vector dx0dt(v0);
        tensor wx
        (
            0,          -omega[2],  omega[1],
            omega[2],   0,          -omega[0],
            -omega[1],  omega[0],   0
        );
        tensor dQdt(wx & Q);

        // dvwdt = body_.calcAcceleration(...)
        scalarField dvwdt(6, 0);
        scalarSquareMatrix Madd(body_.Madd());
        Madd = body_.changeFrame(Madd, Q.T());

        scalarField vw(6); //Matrix for product below to be correct?
        vw[0] = v0[0]; vw[1] = v0[1]; vw[2] = v0[2];
        vw[3] = omega[0]; vw[4] = omega[1]; vw[5] = omega[2];
        scalarSquareMatrix W({6, 6}, 0);
        W.subMatrix(0,0,3,3) = wx; //Note: tensor to Matrix - safe?
        W.subMatrix(3,3,3,3) = wx; //Note: tensor to Matrix - safe?
        scalarField Ftau0 = (Madd * W - W * Madd) * vw;
        vector F0(Ftau0[0], Ftau0[1], Ftau0[2]);
        vector tau0(Ftau0[3], Ftau0[4], Ftau0[5]);

        // Adding -v0 x (A v0 + C omega) term to tau0
        tensor v0x
        (
            0,      -v0[2], v0[1],
            v0[2],  0,      -v0[0],
            -v0[1], v0[0],  0
        );
        scalarSquareMatrix auxMatrix({6, 6}, 0);
        auxMatrix.subMatrix(3,0,3,3) = v0x; //Note: tensor to Matrix - safe?
        scalarField auxField = (auxMatrix * Madd) * vw;
        tau0[0] -= auxField[3];
        tau0[1] -= auxField[4];
        tau0[2] -= auxField[5];

        Info << "F0 = " << F0 << ", tau0 = " << tau0 << endl;
        // Setting body state which is used by calcAcceleration
        body_.setCentreOfRotation(x0);
        body_.setOrientation(Q);
        body_.setVelocity(v0);
        body_.setAngularVelocity(omega);
        dvwdt = body_.calcAcceleration(F0, tau0, DoFs_);
        Info << "dvwdt " << dvwdt << endl;

        for (label n = 0; n<3; n++)
        {
            dYdt[n] = v0[n];
            dYdt[n+12] = dvwdt[n];
            dYdt[n+12+3] = dvwdt[n+3];
        }
        label k = 3;
        for (label m = 0; m<3; m++)
        {
            for (label n = 0; n<3; n++)
            {
                dYdt[k] = dQdt(m,n);
                k++;
            }
        }
    }

    void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
        dfdx = 0.0;
    }
};
