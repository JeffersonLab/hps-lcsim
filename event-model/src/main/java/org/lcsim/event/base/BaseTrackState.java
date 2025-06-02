package org.lcsim.event.base;

import static java.lang.Math.abs;

import java.io.PrintStream;

import org.lcsim.constants.Constants;
import org.lcsim.event.TrackState;

/**
 * Implementation of the org.lcsim.event.TrackState interface.
 * @author Jeremy McCormick
 * @version $Id: BaseTrackState.java,v 1.2 2012/06/18 23:02:14 jeremy Exp $
 */
public class BaseTrackState implements TrackState
{
    public static final int PARAMETERS_SIZE  = 5; // Number of LCIO track parameters.
    public static final int REF_POINT_SIZE   = 3; // Size of reference point (x, y, z).
    public static final int COV_MATRIX_SIZE = 15; // Size of covariance matrix array.
    public static final int MOMENTUM_SIZE   = 3; // Size of momentum (x, y, z).

    private double _def=-666; //use crazy default value
    
    private double[] _parameters = new double[PARAMETERS_SIZE];     // Parameters array.
    private double[] _referencePoint = new double[REF_POINT_SIZE];  // Reference point.
    private double[] _covMatrix = new double[COV_MATRIX_SIZE];      // Covariance matrix.
    private double[] _momentum= new double[MOMENTUM_SIZE];    //momentum array
   
    // local bfield
    private double _blocal=_def;  //initialize with something crazy 
    // Location encoding.
    private int _location = TrackState.AtPerigee; // default location
       
    public BaseTrackState()
    {}
    
    public BaseTrackState(int location)
    {_location = location; }

    
    //fully qualified constructor but with no bfield (so no momentum)
    public BaseTrackState(double[] trackParameters, double[] position, double[] covarianceMatrix, int location)
    {
        _location = location;
        System.arraycopy(trackParameters, 0, _parameters, 0, PARAMETERS_SIZE);
        System.arraycopy(covarianceMatrix, 0, _covMatrix, 0, COV_MATRIX_SIZE);
        System.arraycopy(position,0, _referencePoint, 0, REF_POINT_SIZE);	
    }
    
    // Ctor with parameters and B-field.  
    // The reference point, covariance matrix, and location can be set later.
    public BaseTrackState(double[] parameters, double bfield)
    {
	setBLocal(bfield);
        setParameters(parameters, bfield);
    }
    
    // Fully qualified constructor.
    public BaseTrackState(double[] parameters, double[] referencePoint, double[] covMatrix, int location, double bfield)
    {
	setBLocal(bfield); 
        setParameters(parameters, bfield);
        setReferencePoint(referencePoint);
        setCovMatrix(covMatrix);
        setLocation(location);
	
    }
     
    public int getLocation()
    {
        return _location;
    }
    
    public double[] getReferencePoint()
    {
        return _referencePoint;
    }

    public double[] getCovMatrix()
    {
        return _covMatrix;
    }
    
    public double getD0()
    {
        return _parameters[BaseTrack.D0];
    }

    public double getPhi()
    {
        return _parameters[BaseTrack.PHI];
    }

    public double getZ0()
    {
        return _parameters[BaseTrack.Z0];
    }

    public double getOmega()
    {
        return _parameters[BaseTrack.OMEGA];
    }

    public double getTanLambda()
    {
        return _parameters[BaseTrack.TANLAMBDA];
    }
    public double getBLocal()
    {
        return _blocal;
    }
    
    public void setD0(double d0)
    {
        _parameters[BaseTrack.D0] = d0;
    }

    public void setPhi(double phi)
    {
        _parameters[BaseTrack.PHI] = phi;
    }

    public void setZ0(double z0)
    {
        _parameters[BaseTrack.Z0] = z0;
    }

    public void setOmega(double d)
    {
        _parameters[BaseTrack.OMEGA] = d;
    }

    public void setTanLambda(double d)
    { 
        _parameters[BaseTrack.TANLAMBDA] = d;
    }
    
    public void setLocation(int location)
    {
        if (location > TrackState.LastLocation)
            throw new IllegalArgumentException("The location must be between 0 and " + TrackState.LastLocation);
        this._location = location;
    }

    public void setBLocal(double d)
    { 
        _blocal = d;
    }
    

    // FIXME Should be array copy?
    public void setReferencePoint(double[] referencePoint)
    {
        if (referencePoint.length != REF_POINT_SIZE) throw new IllegalArgumentException("referencePoint.length != " + REF_POINT_SIZE);
        this._referencePoint = referencePoint;
    }
    
    // FIXME Should be array copy?
    public void setCovMatrix(double[] covMatrix)
    {
        if (covMatrix.length != COV_MATRIX_SIZE) throw new IllegalArgumentException("covMatrix.length != " + COV_MATRIX_SIZE);
        this._covMatrix = covMatrix;
    }
    
    // If setParameters or a qualified constructor was not called, this could be null.
    public double[] getMomentum()
    {
	//make sure this has been computed
	if(_momentum[0]== 0 )//it hasn't been computed
	    if(_blocal!=_def)
		computeMomentum();	       	    
        return _momentum;
    }
    
    // Get a track parameter by ordinal.
    public double getParameter(int param)
    {   
        if (param < 0 || param > (PARAMETERS_SIZE - 1))
            throw new IllegalArgumentException("Parameter ordinal " + param + " is invalid.");
        return _parameters[param];
    }
    
    // Get the parameters as a double array.  
    // Use ordinals in BaseTrack for the index into the array.
    public double[] getParameters()
    {
        return _parameters;
    }
    
    /**
     * Set the track parameters.  Computes momentum and charge, also.
     * @param p
     * @param bfield
     */
    public void setParameters(double[] p, double bfield)
    {
	_blocal=bfield;
        copyParameters(p, _parameters);
        computeMomentum();
    }
    
    static final void copyParameters(double[] p1, double[] p2)
    {
        if (p1.length != 5)
            throw new IllegalArgumentException("First array is not size " + PARAMETERS_SIZE);
        if (p2.length != 5)
            throw new IllegalArgumentException("Second aray is not size" + PARAMETERS_SIZE);
        System.arraycopy(p1, 0, p2, 0, PARAMETERS_SIZE);
    }
    
   
    public void computeMomentum(){
	double omega = this.getOmega();
        if(abs(omega) < 0.0000001) omega = 0.0000001;
        double Pt = abs((1./omega) * _blocal * Constants.fieldConversion);  
        double[] momentum = new double[3];
        momentum[0] = Pt * Math.cos(this.getPhi());
        momentum[1] = Pt * Math.sin(this.getPhi());
        momentum[2] = Pt * this.getTanLambda();
	setMomentum(momentum);
    }
    
    public void setMomentum(double[] mom){
	this._momentum[0]=mom[0];
	this._momentum[1]=mom[1];
	this._momentum[2]=mom[2];
    }

     /**
     * Convert object to a String.
     */
    public String toString()
    {
        StringBuffer buff = new StringBuffer();
        buff.append("location = " + getLocation() + "\n");
        buff.append("D0 = " + getD0() + "\n");
        buff.append("phi = " + getPhi() + "\n");
        buff.append("Z0 = " + getZ0() + "\n");
        buff.append("tanLambda = " + getTanLambda() + "\n");
        buff.append("omega = " + getOmega() + "\n");
        buff.append("referencePoint = " + _referencePoint[0] + " " + _referencePoint[1] + " " + _referencePoint[2] + "\n");
	buff.append("bField at ref = "+_blocal+"\n");
        buff.append("covarianceMatrix = ");
        for (int i=0; i<_covMatrix.length; i++) 
        {
            buff.append(_covMatrix[i] + " ");
        }
        buff.append("\n");
        buff.append("momentum = ");
        for (int i=0; i<this.MOMENTUM_SIZE; i++)
        {
            buff.append(_momentum[i] + " ");
        }
        buff.append("\n");
        return buff.toString();
    }
    
    public void printOut(PrintStream ps)
    {
        ps.println(toString());
    }
}
