/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package seqreport;

//ivica res, May05
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import javax.swing.*;
import java.awt.image.*;
import javax.imageio.*;
import javax.swing.ImageIcon;
//import org.jibble.epsgraphics.*;



public class SeqReportEE {

	final static int MARGIN = 35;
	final static int FONTSIZE = 30;
	final static int STRIP_SPACING = 32;
	private float increaseSize = 1;
	private static final double STRIP_WIDTH = 23.158;
//	final float MAXCOLOR = 0.8f;
    final float MAXCOLOR = 1f;  // for Ivana color scheme
	             //float[] x = new float[10];
	ArrayList<Double> coverage;
	ArrayList<Double> altCoverage;
	ArrayList<Integer> resNumbers;
	          // int[] y = new int[10];
	          //    String[] residue = new String[10];
	ArrayList<String> residue;
    float minX = 1000.f;
    float maxX = -1000.f;
  
    private boolean missing = false;

	public SeqReportEE() {
		residue = new ArrayList<String>();
		coverage = new ArrayList<Double>();
		altCoverage = new ArrayList<Double>();
		resNumbers = new ArrayList<Integer>();
	}

	public SeqReportEE(float increaseSize) {
		this.increaseSize = increaseSize;
		residue = new ArrayList<String>();
		coverage = new ArrayList<Double>();
		altCoverage = new ArrayList<Double>();
		resNumbers = new ArrayList<Integer>();
	}
	public void readPoints(File file, int start, int end) {

		int i = 0;
		int cnt = 0;
		int arrayIndex = 0;
		double tempInput = 0;
		double altTempInput = 0;
		int tempExtra = 0;
		String line;
		String stringInput = "";

                String content = null;  
                String ignore = "no";
                String[] too = null;
		try {
                       
			Reader reader = new BufferedReader( new FileReader(file));
                        while ((content = ((BufferedReader)reader).readLine()) != null)  {
                             if ( content.trim().length() == 0 )  {
                              //   System.out.println("zfgafsgasgs");
                                 continue;
                             }
                 //          System.out.println(content);
                            too = content.trim().split("\\s+");
                         //   System.out.println(content);
              //             System.out.println(too[0] + "*" + too[1] + "*" + too[2] + "*" + too[3] + "*");
                            for(String st : too) 
                                
                                if (st.equals(".")) {
                                  ignore = "yes";
           //  System.out.print("*************" + st);

                                }
                          
                           // System.out.println();
                            if (ignore.equals("no")) {
                                coverage.add(Double.parseDouble(too[0]));
                                altCoverage.add(1 - 2 * Double.parseDouble(too[1]));
                                residue.add(too[2]);
                                resNumbers.add(Integer.parseInt(too[3]));
                   //             System.out.println(content);
                            }
                            ignore = "no";
                        }
  
                //       System.exit(0);
                        
                        
                       
			int countMissing = 0;
			for (int j = 0; j < resNumbers.size() - 1; j++) {
				if ( (int)(resNumbers.get(j + 1) - resNumbers.get(j)) != 1) {
					missing = true;
					countMissing ++;
				}
			}
			//	System.out.println("d = " + countMissing);
                        missing =  false;
			if (missing) {
				ArrayList<Double> modCoverage = new ArrayList<Double>();
				ArrayList<Double> modAltCoverage = new ArrayList<Double>();
				ArrayList<Integer> modResnumbers = new ArrayList<Integer>();

				ArrayList<String> modResidue = new ArrayList<String>();
				int k = 0;
				int j = 0;
				System.out.println ("size = " + resNumbers.size());
				for (j = 0; j < resNumbers.size() - 1; j++) {
					modCoverage.add(coverage.get(j));
					modAltCoverage.add(altCoverage.get(j));
					modResidue.add(residue.get(j));
					modResnumbers.add(resNumbers.get(j));


					for (int iii = 0; iii < 3; iii++) {

						if ( (resNumbers.get(j + 1) - resNumbers.get(j)) != 1) {

							modCoverage.add(-10.);
							modAltCoverage.add(-10.);
							modResidue.add("dots");
							modResnumbers.add(-10);
						}
					}

				}
				modCoverage.add(coverage.get(j));
				modAltCoverage.add(altCoverage.get(j));
				modResidue.add(residue.get(j));
				modResnumbers.add(resNumbers.get(j));

				coverage = modCoverage;

				altCoverage = modAltCoverage;
				residue = modResidue;
				resNumbers = modResnumbers;
			}

		}
		catch (Exception e) {
			System.out.println("Exception " + e.getMessage() + " in readRanks()");
		}

	}


	public static void skipLeadingNewLines(StreamTokenizer t) throws IOException {
		while (t.nextToken() == t.TT_EOL);
		t.pushBack();
    }

	class MapPanel extends JComponent {
	//this class does all the drawing

	  	float width = 500;//0;  // frame size
		float height = 500;//0;
		float margin =25;//20;
		float marginTop = 25;//40;
		float squareX = 20;  // legend bar width

		String title;
//		Font numberFont = new Font("Helvetica",Font.PLAIN,20);
//		Font residueFont = new Font("Helvetica",Font.PLAIN,10);
		Font numberFont = new Font("Monospaced",Font.PLAIN,20);
		Font residueFont = new Font("Monospaced",Font.PLAIN,10);


 		MapPanel(String title) {
 			this.title = title;
			setPreferredSize(new Dimension((int)width, (int)height));
 		}

 		MapPanel(String title, float width, float height) {
 			this.title = title;
 			this.width = width;
 			this.height = height;
			setPreferredSize(new Dimension((int)width, (int)height));
 		}

		public Color setPositiveColor(double cvg) {
			float hue = 0.082f;
			float brightness = 1;
			return Color.getHSBColor(hue, (float)cvg, brightness);

		}

		public Color setNegativeColor(double cvg) {
			float hue = 0.667f;
			float brightness = 1;
			return Color.getHSBColor(hue, (float)cvg, brightness);

		}
	    public Color setIvanaColor(double cvg) {
 			//color scheme used by Ivana
 			double red, blue, green, ratio, x;
			int  bin_size;
 			int range = 20;
 			int N = 5;
			int color_index;
 			red = green = blue = 254;

 		//	if (cvg <= 0.05 + 0.00001)
 		//		return new Color(254,(int)(0.83*254), (int)(0.17*254));
 		//	else {

			    color_index = (int) (range * cvg);
			    bin_size = (int)Math.round ( (double)(range - 1) / N);

			    if (cvg <= 0.25) {

				ratio =  (double)(bin_size - color_index + 1)/bin_size;
				//System.out.println ( cvg + " " +  ratio);
                                if(ratio > 1)
                                    ratio = 1;
				red = ratio * 254;
				green = blue = 0.;
				//System.out.println (bin_size +"  " +  color_index );
				//System.exit (1);

			    } else {
				ratio = ( color_index  - range / N) / ( (double) (N-1.)/ N * range);
				if(ratio > 1)
                                    ratio = 1;
                                red = ratio * 254;
				green = blue = red;
			    }
		//	}
			//int re, bl, gr;
			//re = (int)red;
			//bl = (int)blue;
			//gr = (int)green;
			//System.out.println("cvg = " + cvg + " " + re + " " + bl + " " + gr);
			return new Color((int)red, (int)blue, (int)green);

 		}

 		//coloring of data points
 		public Color myColor(double norm) {
 			if (norm < 0.95)     //to avoid starting end ending in red
 				return Color.getHSBColor((float)norm, 1, 1);
 			else
 				return Color.getHSBColor(0.95f, 1, 1);
 		}

		// to draw the color scheme
		public void drawLegend(Graphics2D g, float x, float y) {
			int moveY = 10;
			int i;
			float hue;
     		Color color;
     		Font defaultFont = g.getFont();
     		FontMetrics fm = g.getFontMetrics(numberFont);
     		String top, bottom;

			g.setFont(numberFont);

			//for (i = 0; i < (int)y; i++) {
	//	      for (i = -(int)y; i < (int)y; i++) {
for (i = 0; i < (int)y; i++) {
     			//draw color bar as a series of tiny rectangles
     			hue = MAXCOLOR *i / y;
		//	   	if (hue < 0.)
			//		color = setNegativeColor(-hue);
			//	   	else
				//		color = setPositiveColor(hue);
				
    					color = setIvanaColor(hue);
// 	 			//	color = myColor(hue);
//				if (i<=0)
				Rectangle2D.Float rectangle;
 				g.setPaint(color);
				if (i < 0)
					rectangle = new Rectangle2D.Float(x, (float)y - i + moveY,squareX,1);
				else
					rectangle = new Rectangle2D.Float(x, (float)y - i + moveY,squareX,1);
 		       	g.fill(rectangle);
   			}
   			   			g.setPaint(Color.black);
 			g.drawString("0",x+squareX+5,i + 6 + moveY);
			g.drawString("-1",x+squareX+5,i * 2 + 6+ moveY);
      		i=0;
      //		g.setPaint(myColor(i));
      //		g.setPaint(setIvanaColor(i));

 			g.drawString("1",x+squareX+5,i + 6 + moveY);

  			g.setPaint(Color.black);
   			g.setFont(defaultFont);

  		}


   		public void paintComponent(Graphics _g) {
   			drawSequence(_g,getWidth(), getHeight(),0);
   		}

   		public void drawSequence(Graphics _g, float width, float height, int shift) {
 		//put all this together and paint
			float proteinWidth = (float)STRIP_WIDTH * increaseSize;
 			int i, j, titleCenter;
 			int position = 0;
 			int max = coverage.size()/50 + 1;
   	   		float delta, stripe, lineWidth;
			float cov=1.1f;//0;
			float spacing;

//			Font titleFont = new Font("Helvetica",Font.PLAIN,30);
			Font titleFont = new Font("Monospaced",Font.PLAIN,30);

			FontMetrics fm = _g.getFontMetrics(titleFont);
		    Rectangle2D.Float rectangle, anotherRectangle;
		    GeneralPath boundary = new GeneralPath();
			Graphics2D g = (Graphics2D)_g;

			g.setColor(Color.black);
			g.setFont(titleFont);
			titleCenter = fm.stringWidth(title);
//			g.drawString(title,(int)(width / 2 - titleCenter/2), (int)(marginTop-15));
		//	stripe = (width - 4 * margin) / (double)coverage.size();
			stripe = (width - 2 * margin) / 50;
			spacing = (height - margin - marginTop) / 19;
			//	STRIP_WIDTH = 23.158f;  //former spacing
			//getting nice lines
   		 	g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
						       RenderingHints.VALUE_ANTIALIAS_ON);
			g.setStroke(new BasicStroke(1.0f));

				delta = 0;int once = 0;
				int countResidues = 0;
    		for (i = 0; i < coverage.size(); i++) {
			   		if (i%50 == 0)
		     			delta= 0;

//		 	 		g.setPaint(myColor( coverage.get(i)));
		 	 		g.setPaint(setIvanaColor(coverage.get(i)));
                                        
              //            g.setPaint(setIvanaColor(0.05));
                           //             System.out.println("cov " + i + " " + coverage.get(i) );
		  			if(coverage.get(i) / MAXCOLOR > cov)
		  	  			g.setPaint(Color.blue);
			      	if (i != coverage.size()-1 ||  ((i+1)%50 == 0) )
					    lineWidth = stripe*1.1f;//1.25f; //stripe width slightly enlarged
		  	  								  //to avoid white stripes between adjacent stripes
		  	  		else
		  	  			lineWidth = stripe;
		  	  		if ( (i != coverage.size()-1) && coverage.get(i + 1) < 0)
		  	  			lineWidth = stripe;

					rectangle = new Rectangle2D.Float(margin + delta, marginTop + position,lineWidth,
													  proteinWidth/2);
					anotherRectangle = new Rectangle2D.Float(margin + delta, marginTop + position +
															proteinWidth/2 ,lineWidth,proteinWidth/2);

					if(coverage.get(i) >= 0) {
   		       			g.fill(rectangle);
						//	g.setPaint(setIvanaColor(altCoverage.get(i)));
						g.setPaint(setNegativeColor(altCoverage.get(i)));
						if (altCoverage.get(i) < 0)
							g.setPaint(setNegativeColor(-altCoverage.get(i)));
						else
							g.setPaint(setPositiveColor(altCoverage.get(i)));
						g.fill(anotherRectangle);
 	       				g.setPaint(Color.black);
						boundary.moveTo(margin + delta, marginTop + position);
						boundary.lineTo(margin + delta + lineWidth, marginTop + position);
	       				boundary.moveTo(margin + delta, marginTop + position + proteinWidth/2);
	       				boundary.lineTo(margin + delta + lineWidth, marginTop + position + proteinWidth/2);
	       				boundary.moveTo(margin + delta, marginTop + position + proteinWidth);
	       				boundary.lineTo(margin + delta + lineWidth, marginTop + position + proteinWidth);



						g.setPaint(Color.black);

				      	g.setFont(residueFont);
						String axes = "";
						axes = residue.get(i);
				//	if (i < x.length - 5)
						//	g.drawString(axes,margin + delta + 1.5f, marginTop + position + stripe + 25) ;
						g.drawString(axes,margin + delta + 1.5f, marginTop + position + proteinWidth+
									 (int)(FONTSIZE*.4)) ;
			      		if(countResidues%10 == 0){// || (countResidues+1)%10 == 0) {
//			      		if (i%50== 0 || (i+1)%10 == 0) {
				     		axes = Integer.toString(resNumbers.get(i));//Integer.toString(shift + i+1);
				     		g.drawString(axes,margin + delta, marginTop + position +
										 proteinWidth + 2*(int)(FONTSIZE*.4));

				     	}
				     	countResidues++;
				    }
				    else {
				    	Ellipse2D circle = new Ellipse2D.Float(margin + delta + lineWidth / 2 - 1.5f,
															   marginTop + position + proteinWidth/2, 3, 3);
				    	g.setPaint(Color.black);
				    	g.fill(circle);


				    }

 			  		delta +=stripe;



 					if ((i+1)%50 == 0) {

 						position += proteinWidth + STRIP_SPACING; // seq width + space between sequences
 					}

 				}



				g.draw(boundary);
				g.translate(width - 2*margin, marginTop);

    		}
		}   //end of class MapPanel



    public static void main(String[] args) throws Exception {
		if (args.length < 2) {
                    
                     System.out.println("Usage: java SeqReportEE <coverage file> \"figure name\" ");
                   //  System.out.println("Usage: java SeqReportEE <coverage file> \"figure name\" [scaling factor(default 1)]");
 			System.out.print("<coverage file> is a 4 column file ");
 			System.out.println("of the form \"coverage_conservation coverage_specialization residue_name residue_number\" ");
// 			System.out.println("If no numbers given, the whole coverage file will be used");
// 			System.out.println("If <end number> not given, EOF is assumed to be the end");
 			System.exit(0);
		}

		int w, h, shift, end;
		float dimX, dimY;
    	Graphics2D g;
        BufferedImage image;

		// create a new Plot object
		SeqReportEE plot;// = new SeqReportEE();

		// following options removed
		// So now only 3 args max
		// the third one increases the strip width
		shift = 0;
		end = 10000;

// 		//read input file (should have 3 columns: x,y and z
// 		if (args.length >= 3 && Integer.parseInt(args[2]) > 0) {
// 			shift = Integer.parseInt(args[2]) - 1;
// 		}
// 		else
// 		    shift = 0;
// 		if (args.length == 4 && Integer.parseInt(args[3]) > Integer.parseInt(args[2]))
// 		    end = Integer.parseInt(args[3]) - 1;
// 		else
// 		    end = 10000;
		if (args.length == 3) {
			plot = new SeqReportEE(Float.parseFloat(args[2]));
		}
		else
			plot = new SeqReportEE();

       plot.readPoints(new File(args[0]), shift, end);
		//panel used for drawing
     //   MapPanel map = plot.new MapPanel(args[1]);


	int nmbOfStripes = plot.coverage.size()/ 50 + 1;
//	System.out.println(plot.coverage.size());

dimY = (float)(nmbOfStripes * STRIP_WIDTH * plot.increaseSize +  (nmbOfStripes - 1)* STRIP_SPACING) + 2*MARGIN;

	//  dimY *= 23.157f;
	//	  dimY += 60;
	//dimY += 76;
	//dimY += 50;

  //   dimY = 500;

//	System.out.println(dimY);
	MapPanel map = plot.new MapPanel(args[1], 500, dimY);

        w = (int)map.width;
        h = (int)map.height;

		// produce a png image
       	image = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
        g = image.createGraphics();
       	g.setColor(Color.white);
       	g.fillRect(0, 0, w, h);
		map.drawSequence(g,w,h,shift);
		//	map.drawLegend(g,10,100);
       	ImageIO.write(image, "png", new File(args[1]+".png"));

		//produce an eps image
		//       g = new EpsGraphics2D("map");
		//    map.drawSequence(g,w,h,shift);
        //PrintWriter out = new PrintWriter(new FileWriter(args[1]+".eps"));
        //out.println(g.toString());
        //out.close();

       //create MapPanel JComponent and give it to TreeFrame
	//	new TreeFrame(plot.new MapPanel(args[1],500, dimY));

   }

}   //end of class Plot

class TreeFrame extends JFrame {
	// JFrame which draws JComponent
    TreeFrame(JComponent cmp) {
    	getContentPane().setBackground(Color.white);
   		getContentPane().add(cmp);
   	    pack();
   	    addWindowListener(new WindowAdapter(){
    		public void windowClosing(WindowEvent evt){
				System.exit(0);
      		}
    	});
    	setVisible(true);
    }
}
