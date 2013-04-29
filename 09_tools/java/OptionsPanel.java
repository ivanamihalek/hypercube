/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package seqreport;

/**
 *
 * @author Ivica 2009
 */

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import static seqreport.GUILabels.*;
import static java.awt.GridBagConstraints.*;

public class OptionsPanel extends JPanel implements ActionListener {
    private JButton snp, fragment, gene, runButton;
    private JLabel label;
    private JTextField field, scaling;
    private static final int COLUMNS = 10;
    private File snpInfo, geneFragment, wholeGene;
    private String outputFileName = OUTPUT_FILE_NAME;
    private String scalingFactor = "1";
    private String[] options = new String[3];

    public static void main(String[] args) {

        new OptionsPanel();
    }

    public static void show(JPanel panel) {
        JFrame frame = new JFrame();
        frame.setSize(300, 200);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(panel);
        frame.setVisible(true);
    }

    public OptionsPanel() {
        setName(NAME);
        label = new JLabel(OUTPUT_LABEL_TEXT);
        initialize();
        show(this);
    }

    private void initialize() {
        snp = new JButton(SNP_BUTTON);
        snp.addActionListener(this);
        fragment = new JButton(FRAGMENT_BUTTON);
        fragment.addActionListener(this);
        gene = new JButton(GENE_BUTTON);
        gene.addActionListener(this);
        runButton = new JButton(RUN_BUTTON);
        runButton.addActionListener(this);
        field = new JTextField(OUTPUT_FILE_NAME, COLUMNS);
        field.addActionListener(this);
        scaling = new JTextField("1");
        scaling.addActionListener(this);
        createLayout();
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == snp) {
            snpInfo = loadFile();
            options[0] = snpInfo.toString();
 
            
        }
        if (e.getSource() == scaling) {
            scalingFactor = scaling.getText().trim();
            options[2] = scalingFactor;
        }
//            geneFragment = loadFile();
//        }
//        if (e.getSource() == gene) {
//            wholeGene = loadFile();
//        }
        if (e.getSource() == field) {
            String output =  field.getText().trim();
            options[1] = output;
//            if (output.isEmpty() == false) {
//                outputFileName = output;
//            }
//            else
//                outputFileName = OUTPUT_FILE_NAME;
        }
        if (e.getSource() == runButton) {
            //start the application if all required files are present
//            Object[] objects = {snpInfo, geneFragment};
      //      if (snpInfo != null) {
                try{
String output =  field.getText().trim();
System.out.println(output);options[1] = output;options[2] = scalingFactor;
                    SeqReportEE.main(options);
                }
                catch(Exception ex) {
                    System.out.println("ex = " + ex.getMessage());
                }
//            }
            }}//}


    private File loadFile() {
        JFileChooser chooser = new JFileChooser( );
        chooser.showOpenDialog(this);
        return chooser.getSelectedFile();
    }

    private boolean start(Object[] objects) {
        for (int i = 0; i < objects.length; i++) {
            if (objects[i] == null) {
                System.out.println("i = " + i);
                return false;
            }
        }
        return true;
    }

    // the rest of the code is the look and feel of OptionsPanel
    private void createLayout() {
        label = new JLabel(OUTPUT_LABEL_TEXT);
        setLayout(new BorderLayout());
        JPanel panel = createCentralPanel(snp, fragment,gene, runButton);
        add(panel, BorderLayout.CENTER);
    }

    private JPanel createCentralPanel(JButton snp, JButton fragment, JButton gene,
                             JButton run ) {
        JPanel panel = new JPanel();
        panel.setLayout(new BoxLayout(panel, BoxLayout.PAGE_AXIS));

        panel.add(Box.createRigidArea(new Dimension(0, 6)));
        snp.setAlignmentX(Component.LEFT_ALIGNMENT);
        fragment.setAlignmentX(Component.LEFT_ALIGNMENT);
        gene.setAlignmentX(Component.LEFT_ALIGNMENT);
        panel.add(snp, BorderLayout.NORTH);
        panel.add(fragment, BorderLayout.NORTH);
//        panel.add(gene, BorderLayout.NORTH);
    
//        panel.add(Box.createRigidArea(new Dimension(0, 6)));
        
        panel.add(createFieldsPanel(label, field), BorderLayout.SOUTH);
 
        panel.add(new JLabel("Scale Factor"), BorderLayout.NORTH);
 //       panel.add(new JTextField("1"), BorderLayout.NORTH);
        panel.add(scaling);
               panel.add(run, BorderLayout.SOUTH);
        panel.setBorder(BorderFactory.createEmptyBorder(8, 8, 8, 8));
        return panel;
    }

    private JPanel createFieldsPanel(JLabel label, JTextField field) {
        GridBagLayout layout = new GridBagLayout();
        JPanel panel = new JPanel(layout);
        addField(panel, layout, 0, label, field);
            return panel;
    }

    private void addField(JPanel panel, GridBagLayout layout, int row,
                        JLabel label, JTextField field) {
        Insets insets = new Insets(3, 3, 3, 3);
        layout.setConstraints(label,
            new GridBagConstraints(
                0, row,  // x, y
                1, 1,  // gridwidth, gridheight
                40, 1, // weightx, weighty
                LINE_END,
                NONE, // fill
                insets, 0, 0)); // padx, ipady
        layout.setConstraints(field,
            new GridBagConstraints(1, row,
                2, 1, 60, 1, CENTER, HORIZONTAL,
                insets, 0, 0));
        layout.setConstraints(label,
            new GridBagConstraints(0, row, 1, 1, 40, 1,
                LINE_END, NONE, insets, 0, 0));
        panel.add(label);
        panel.add(field);
   }
}