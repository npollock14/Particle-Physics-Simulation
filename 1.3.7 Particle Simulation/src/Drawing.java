import java.awt.Canvas;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;

public class Drawing extends Canvas implements MouseListener, MouseMotionListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = -5082417085041890736L;
	ArrayList<Double> data = new ArrayList<Double>();
	double xLow, xHi, yLow, yHi;
	int width, height;
	int radius = 2;
	int radiusMult = 1;
	double selectedY;

	public Drawing(int width, int height, ArrayList<Double> data) {
		this.data = clean(data);
		this.width = width;
		this.height = height;
		this.addMouseListener(this);
		this.addMouseMotionListener(this);

		this.radius = (int) Math.sqrt((radiusMult * 100000 / data.size()));
		if (this.radius < 2) {
			this.radius = 2;
		}
		if(this.radius > 10) {
			this.radius = 10;
		}

	}

	private ArrayList<Double> clean(ArrayList<Double> foo) {
		for (int i = 0; i < foo.size(); i++) {
			if (foo.get(i).isNaN() || foo.get(i) == null) {
				foo.remove(i);
			}
		}
		return foo;
	}

	public void paint(Graphics g) {
		g.setFont(new Font("TimesRoman", Font.BOLD, 15)); 
		yLow = data.get(0);
		for (int i = 0; i < data.size(); i++) {
			if (data.get(i) < yLow) {
				yLow = data.get(i);
			}
		}
		yHi = data.get(0);
		for (int k = 0; k < data.size(); k++) {
			if (data.get(k) > yHi) {
				yHi = data.get(k);
			}
		}
		xLow = 0;
		xHi = data.size() - 1;
		
		if (yHi > 0 && yLow < 0) {
			g.drawLine(0, (int) (height - (-yLow * (height / (yHi - yLow)))), width,
					(int) (height - (-yLow * (height / (yHi - yLow)))));
		}

		g.drawString("High: " + yHi, 0, 15);
		g.drawString("Low: " + yLow, 0, 30);
		
		for (int l = 0; l < data.size(); l++) {
			g.fillOval((int) (l * width / xHi), (int) (height - ((data.get(l) - yLow) * (height / (yHi - yLow)))), 1,
					1);
		}
		
		

	}
	
	public void mousePressed(MouseEvent e) {
		Graphics g = this.getGraphics();
		this.getGraphics().clearRect(0, 0, width, height);
		paint(this.getGraphics());
		this.getGraphics().drawLine(e.getX(), height, e.getX(), 0);
		g.setColor(Color.RED);
		g.setFont(new Font("TimesRoman", Font.BOLD, 15)); 
		g.drawString("Selection: " + getY(e.getX()), 0, 45);
		g.setColor(Color.black);
	}

	public void mouseReleased(MouseEvent e) {
		this.getGraphics().clearRect(0, 0, width, height);
		paint(this.getGraphics());
	}

	public void mouseEntered(MouseEvent e) {
	   
	}

	public void mouseExited(MouseEvent e) {
		
	}

	public void mouseClicked(MouseEvent e) {
		
	}

	private String getY(int Px) {
		int x = (int)(Px*data.size()/width);
		
		return data.get(x).toString();
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		Graphics g = this.getGraphics();
		this.getGraphics().clearRect(0, 0, width, height);
		paint(this.getGraphics());
		this.getGraphics().drawLine(e.getX(), height, e.getX(), 0);
		g.setColor(Color.RED);
		g.setFont(new Font("TimesRoman", Font.BOLD, 15)); 
		g.drawString("Selection: " + getY(e.getX()), 0, 45);
		g.setColor(Color.black);
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		
		
	}

}