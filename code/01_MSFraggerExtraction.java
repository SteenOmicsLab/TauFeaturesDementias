import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class MSFraggerExtraction {
	
	private static String convertmodsfrommodseq(String modseq) {
		modseq = modseq.replaceAll("n\\[42\\.0106\\]","");
		int idx = 0;
		boolean inmod = false;
		String modstring = "";
		for(int i = 0; i < modseq.length(); ++i) {
			if(inmod==false && modseq.charAt(i)=='[') {
				inmod = true;
				if(!modstring.equals("")) {
					modstring = modstring + ",";
				}
				modstring = modstring + idx + modseq.charAt(i-1) + "(";
			} else if(inmod==true && modseq.charAt(i)!=']') {
				modstring = modstring + modseq.charAt(i);
			} else if(inmod==true && modseq.charAt(i)==']') {
				modstring = modstring + ")";
				inmod = false;
			} else {
				idx = idx + 1;
			}
		}
		return modstring;
	}
	
	private static String reformmods(String modifications, Map<String,String> modmap) {
		String modstring = "";
		modifications = convertmodsfrommodseq(modifications);
		if(!modifications.equals("")) {
			String[] allmodifications = modifications.split(",");
			for(int i = 0; i<allmodifications.length; ++i) {
				String mod = allmodifications[i];
				mod = mod.replaceAll("\\)","");
				String[] modsplit = mod.split("\\(");
				String aa = modsplit[0].substring(modsplit[0].length()-1);
				String site = modsplit[0];
				site = site.replaceAll(aa,"");
				String modification = modmap.get(modsplit[1]);
				if(!modstring.equals("")) {
					modstring = modstring + "; ";
				}
				modstring = modstring + modification + "(" + aa + ")@" + site;
			}
		}
		return modstring;
	}
	
	private static String correctingMods(String mods, int pepstart, String pep, String accession) {
		String newmodslist = "";
		String[] protmodlist = mods.split(";");
		for(int i = 0; i < protmodlist.length; ++i) {
			String mod = protmodlist[i].trim();
			int idx = mod.indexOf("@");
			String mod_name = "";
			if(idx>-1) {
				mod_name= mod.substring(0, idx);
				int idx2 = mod_name.indexOf(")");
				if(idx2>-1) {
					mod_name = mod_name.substring(0, idx2+1);
				}
				String siteidx = mod.substring(idx+1);
				int siteidx_i = -1;
				if(siteidx.startsWith("N-term")) {
					siteidx_i = 0;
				} else if(siteidx.startsWith("C-term")){
					siteidx_i = pep.length();
				} else {
					siteidx_i = Integer.parseInt(siteidx);
				}
				siteidx_i = siteidx_i + pepstart;
				
				if(mod_name.equals("Phospho(S)") || mod_name.equals("Phospho(T)") || mod_name.equals("Phospho(Y)") || 
						/*mod_name.equals("Acetyl") || */mod_name.equals("Acetyl(K)") || mod_name.equals("GG(K)") || mod_name.equals("GG(S)") ||
						mod_name.equals("Methyl(K)") || mod_name.equals("Methyl(R)") || mod_name.equals("Oxidation(M)") || 
						mod_name.equals("Deamidated(N)") || mod_name.equals("Deamidated(Q)") || mod_name.equals("Oxidation(P)") ||
						mod_name.equals("Deamidated(R)") || mod_name.equals("Citrullination(R)")) {
					if(mod_name.equals("Deamidated(R)")) {
						mod_name = "Citrullination(R)";
					}
					if(accession.equals("3R")) {
						if(siteidx_i>20) {
							siteidx_i = siteidx_i + 31;
						}
						siteidx_i = siteidx_i + 254;
					} else if(accession.equals("0N")) {
						if(siteidx_i>20) {
							siteidx_i = siteidx_i + 58;
						}
						siteidx_i = siteidx_i + 24;
					} else if(accession.equals("1N")) {
						if(siteidx_i>49) {
							siteidx_i = siteidx_i + 29;
						}
						siteidx_i = siteidx_i + 24;
					} else if(accession.equals("epsilon")) {
						siteidx_i = siteidx_i + 376;
					} else if(accession.equals("Ex3")) {
						if(siteidx_i>7) {
							siteidx_i = siteidx_i + 14;
						}
						siteidx_i = siteidx_i + 32;
					} else if(accession.equals("Ex5")) {
						if(siteidx_i>7) {
							siteidx_i = siteidx_i + 28;
						}
						siteidx_i = siteidx_i + 96;
					}
					String tmp = mod_name + "@" + siteidx_i;
					if(pepstart>-1) {
						if(!newmodslist.equals("")) {
							newmodslist = newmodslist + "; ";
						}
						newmodslist = newmodslist + tmp;
					}
				}
			}
		}
		return newmodslist;
	}
	
	private static Integer offsetfromAccession(String accession, int siteidx_i) {
		int offset = 0;
		if(accession.equals("3R")) {
			if(siteidx_i>20) {
				offset = offset + 31;
			}
			offset = offset + 254;
		} else if(accession.equals("0N")) {
			if(siteidx_i>20) {
				offset = offset + 58;
			}
			offset = offset + 24;
		} else if(accession.equals("1N")) {
			if(siteidx_i>49) {
				offset = offset + 29;
			}
			offset = offset + 24;
		} else if(accession.equals("epsilon")) {
			offset = offset + 376;
		} else if(accession.equals("Ex3")) {
			if(siteidx_i>7) {
				offset = offset + 14;
			}
			offset = offset + 32;
		} else if(accession.equals("Ex5")) {
			if(siteidx_i>7) {
				offset = offset + 28;
			}
			offset = offset + 96;
		}
		return offset;
	}
	
	
	private static String checkNtermCleavage(String sequence, String peptide, String accession) {
		boolean cleavage = false;
		String cleavagestring = "";
		int index = sequence.indexOf(peptide);
		String AA_first = peptide.substring(0,1);
		String AA_before = "";
		if(index>0) {
			AA_before = sequence.substring(index-1,index);
		}
		if(!AA_before.equals("K") && !AA_before.equals("R")) {
			cleavage = true;
		}
		if(cleavage) {
			cleavagestring = "N-term(" + AA_first + ")@" + (index + 1 + offsetfromAccession(accession,index+1));
		}
		return cleavagestring;
	}
	
	private static String checkCtermCleavage(String sequence, String peptide, String accession) {
		boolean cleavage = false;
		int index = sequence.indexOf(peptide) + peptide.length();
		String cleavagestring = "";
		String AA_last = peptide.substring(peptide.length()-1);
		if(!AA_last.equals("K") && !AA_last.equals("R")) {
			cleavage = true;
		}
		if(cleavage) {
			cleavagestring = "C-term(" + AA_last + ")@" + (index + offsetfromAccession(accession,index));
		}
		return cleavagestring;
	}
	
	private static boolean checkmods(String ProtModifications) {
		boolean result = false;
		
		Set<String> protmodtable_mods = new HashSet<String>();
		
		String[] protmodlist = ProtModifications.split(";");
		for(int i = 0; i < protmodlist.length; ++i) {
			String mod = protmodlist[i].trim();
			int idx = mod.indexOf("@");
			String mod_name = "";
			if(idx>-1) {
				mod_name= mod.substring(0, idx);
				int idx2 = mod_name.indexOf(")");
				if(idx2>-1) {
					mod_name = mod_name.substring(0, idx2+1);
				}
				if(mod_name.equals("Deamidated(R)")) {
					mod_name = "Citrullination(R)";
				}
				protmodtable_mods.add(mod_name);
			}
		}
		int count = 0;
		int count_interest = 0;
		if(protmodtable_mods.contains("Phospho(S)") || protmodtable_mods.contains("Phospho(T)") || protmodtable_mods.contains("Phospho(Y)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("GG(K)") || protmodtable_mods.contains("GG(S)")) {
			count += 1;
			count_interest += 1;
		}
		if (/*protmodtable_mods.contains("Acetyl") || */protmodtable_mods.contains("Acetyl(K)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Methyl(K)") | protmodtable_mods.contains("Methyl(R)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Citrullination(R)") | protmodtable_mods.contains("Deamidated(R)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Deamidated(N)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Oxidation(P)")) {
			count += 1;
			count_interest += 1;
		}
		if (protmodtable_mods.contains("Oxidation(M)")) {
			count += 1;
		}
		if (protmodtable_mods.contains("Deamidated(Q)")) {
			count += 1;
		}
		
		if(count==protmodtable_mods.size() && count_interest>0) {
			result = true;
		} else {
			result = false;
		}
		
		return result;
	}
	
	private static String addPTMsequences(String sequence, String modifications, String ProtModifications) {
		List<String> protmodtable_mods = new ArrayList<String>();
		List<String> protmodtable_name = new ArrayList<String>();
		List<Integer> protmodtable_idxs = new ArrayList<Integer>();
		String modsequence = "";
//		System.out.println(sequence + "\t" + modifications + "\t" + ProtModifications);
		if(!sequence.equals("") && !modifications.equals("") && !ProtModifications.equals("")) {
			modsequence = sequence;
			String[] modlist = modifications.split(";");
			String[] protmodlist = ProtModifications.split(";");
			for(int j = protmodlist.length-1; j >= 0; --j) {
				String mod = protmodlist[j].trim();
				if(!mod.contains("-term")) {
					int idx = mod.indexOf("(");
					String mod_name = "";
					if(idx>-1) {
						mod_name= mod.substring(0, idx);
						int idx2 = mod.indexOf(")");
						String mod_site = mod.substring(idx+1, idx2);
						int idx3 = mod.indexOf("@");
						int mod_index = Integer.parseInt(mod.substring(idx3+1)) - 1;
						protmodtable_mods.add(mod_name + "(" + mod_site + ")");
						protmodtable_name.add(mod_name);
						protmodtable_idxs.add(mod_index);
					}
					
				} else {
					int idx = mod.indexOf("@");
					String mod_name = mod.substring(0, idx);
					String mod_site = mod.substring(idx+1);
					if(mod_site.equals("N-term")) {
						protmodtable_mods.add(mod_name + "(" + mod_site + ")");
						protmodtable_name.add(mod_name);
						protmodtable_idxs.add(-1);
					} else if(mod_site.equals("C-term")) {
						protmodtable_mods.add(mod_name + "(" + mod_site + ")");
						protmodtable_name.add(mod_name);
						protmodtable_idxs.add(sequence.length()+1);
					}
				}
			}			
			
			for(int i = modlist.length-1; i >= 0; --i) {
				String mod = modlist[i].trim();
				if(!mod.contains("-term")) {
					int idx = mod.indexOf("(");
					String mod_name = mod.substring(0, idx);
					int idx2 = mod.indexOf(")");
					String mod_site = mod.substring(idx+1, idx2);
					int idx3 = mod.indexOf("@");
					int mod_index = Integer.parseInt(mod.substring(idx3+1)) - 1;
					boolean found = false;
					int protmodidx = 0;
					while(found!=true && protmodidx<protmodtable_mods.size()) {
						if(protmodtable_name.get(protmodidx).equals(mod_name)) {
							if(protmodtable_mods.get(protmodidx).equals(mod_name + "(" + mod_site + ")")) {
								found = true;
							} else if(mod_site.equals("N-term") || mod_site.equals("C-term")){
								found = true;
							}
						}
						protmodidx += 1;
					}
					if(found==true) {
						protmodidx = protmodidx - 1;
						System.out.println(modsequence);
						System.out.println(mod_index);
						if(modsequence.equals("HVPGGGSVQIVYKPVDLSKVTSK")) {
							boolean test = true;
						}
						String pepsite = "" + modsequence.charAt(mod_index);
						if(mod_site.equals(pepsite)) {
							modsequence = modsequence.substring(0, mod_index + 1) + "(" + mod_name + ")" + modsequence.substring(mod_index + 1);
//							String protmod_mod = protmodtable_mods.remove(protmodidx);
//							String protmod_nam = protmodtable_name.remove(protmodidx);
//							Integer protmod_idx = protmodtable_idxs.remove(protmodidx);
						}
					}
				} else {
					int idx = mod.indexOf("@");
					String mod_name = mod.substring(0, idx);
					String mod_site = mod.substring(idx+1);
					boolean found = false;
					int protmodidx = 0;
					while(found!=true && protmodidx<protmodtable_mods.size()) {
						if(protmodtable_name.get(protmodidx).equals(mod_name)) {
							if(protmodtable_mods.get(protmodidx).equals(mod_name + "(" + mod_site + ")")) {
								found = true;
							} else if(mod_site.equals("N-term") || mod_site.equals("C-term")){
								found = true;
							}
						}
						protmodidx += 1;
					}
					if(found==true) {
						protmodidx = protmodidx - 1;
						if(mod_site.equals("N-term")) {
							modsequence = "(" + mod_name + ")" + modsequence;
//							String protmod_mod = protmodtable_mods.remove(protmodidx);
//							String protmod_nam = protmodtable_name.remove(protmodidx);
//							Integer protmod_idx = protmodtable_idxs.remove(protmodidx);
						} else if(mod_site.equals("C-term")) {
							modsequence = modsequence + "(" + mod_name + ")";
//							String protmod_mod = protmodtable_mods.remove(protmodidx);
//							String protmod_nam = protmodtable_name.remove(protmodidx);
//							Integer protmod_idx = protmodtable_idxs.remove(protmodidx);
						}
					}
				}
			}
		} else if(!sequence.equals("") && modifications.equals("")) {
			modsequence = sequence;
		}
		return modsequence;
	}
	
	

	public static void main(String[] args) {
		String taufile = "D:/Work/Tau/uniprot_tau_isoforms2N4R-0N-1N-3R.fasta";//args[0]; // "D:/Data/GFAP/AlignmentData/AlignmentDBs/uniprot_tau_isoforms2N4R-0N-1N-3R.fasta";
    	List<String> headers_ordered = new ArrayList<String>();
    	Map<String, String> taumap = new HashMap<String,String>();
    	
    	BufferedReader reader;
    	try {
    		reader = new BufferedReader(new FileReader(new File(taufile)));
    		String line = "";
    		String header = "";
    		String sequence = "";
    		while((line=reader.readLine())!=null) {
    			if(line.startsWith(">")) {
    				if(!header.equals("")) {
    					headers_ordered.add(header);
    					taumap.put(header, sequence);
    					header = "";
    					sequence = "";
    				}
    				header = line.split("\\|")[1];
    			} else {
    				sequence = sequence + line.trim();
    			}
    		}
    		headers_ordered.add(header);
			taumap.put(header, sequence);
    		reader.close();
    	} catch (IOException e) {
    		// TODO auto-generated catch block
    		e.printStackTrace();
    	}
    	
//    	System.out.println(headers_ordered);
//    	System.out.println(taumap);
//    	System.out.println(taumap.get(headers_ordered.get(0)) + "\n" + "AKTDHGAEIVYK");
//    	System.out.println(taumap.get(headers_ordered.get(0)).contains("AKTDHGAEIVYK"));
//    	System.out.println(taumap.get(headers_ordered.get(0)).indexOf("AKTDHGAEIVYK"));
    	
    	
        String filetablelocation = "D:/Work/Tau/00_FDR_files_MSFragger_WH.txt";//args[1]; // "D:/Data/Simon_HMW_paper/Addition_Mascot/00_FDR_files_LMW.txt";
    	String path = "D:/Work/Tau/";//args[2]; // "D:/Data/Simon_HMW_paper/Addition_Mascot/";
    	String outputfile = "D:/Work/Tau/MSFragger_WH_Insol2_wInt_summary_Int.csv";//args[3]; // "D:/Data/Simon_HMW_paper/Addition_Mascot/MASCOT-LMW_cutoff_0-05FDR_TAU.csv";
//    	double threshold = 0.01;//Double.parseDouble(args[4]); //1; //
    	String id = "P10636";//args[4]; // "P10636"; 
//    	int fdrtype = (int)Double.parseDouble(args[6]); // 17 == global FDR; 18 == local FDR
    	String sep = ",";
    	List<String> files = new ArrayList<String>();
    	
    	try {
			reader = new BufferedReader(new FileReader(new File(filetablelocation)));
			String line = "";
			while((line=reader.readLine())!=null) {
				files.add(line);
			}
			reader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
    	
    	
    	Writer writer;
		try {
			writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputfile)));
//			writer.write("Groupfile" + sep + "FDR" + sep + "FDR Threshold" + sep + "Conf Threshold\n");
			writer.write("Sample" + sep + "Accessions" + sep + "Conf" + sep + "Sequence" + sep + "Modificaitons" + sep + "CorrectedModifications" + sep + "ProteinModifications" + sep + "ModifiedSequence" + sep + "Cleavages" + sep + "CleavageSites" + sep + "dMass" + sep + "Obs MW" + sep + "Obs m/z" + sep + "Theor MW" + sep + "Theor m/z" + sep + "Theor z" + sep + "Intensity" + "\n");
			
			for(String file : files) {
		    	
		    	String fileLocation = path + File.separatorChar + file;
//		    	String[] filesplit = file.split("/");
//		    	String newfile = filesplit[filesplit.length-1];
//		    	String[] newfilesplit = newfile.split("_");
//		    	String newfilename = newfilesplit[0] + "_" + newfilesplit[2];
		    	
//		    	String group = newfilesplit[0];
//		    	String region = newfilesplit[1];
//		    	String sample = newfilesplit[2];
//		    	String fraction = newfilesplit[3];
//		    	if(group.equals("CTRL")) {
//		    		group = "Control";
//		    	}
//		    	if(region.equals("PF")) {
//		    		region = "Prefrontal";
//		    	}
		    	
//		    	String newfilename = newfile;
//		    	String newfilename = newfilesplit[2] + "_" + newfilesplit[4];
		    	
		    	Map<String,String> modlookup = new HashMap<String,String>();
		    	modlookup.put("15.9949","Oxidation");
		    	modlookup.put("42.0106","Acetyl");
		    	modlookup.put("71.0371","Propionamide");
		    	modlookup.put("79.9663","Phospho");
		    	modlookup.put("114.0429","GG");
		    	int prot_acc_idx = -1;
		    	int pep_seq_idx = -1;
		    	int pep_mod_seq_idx = -1;
		    	int pep_exp_z_idx = -1;
		    	List<Integer> exp_idxs = new ArrayList<Integer>();
		    	List<String> exp_names = new ArrayList<String>();
		    	
				BufferedReader pepsumreader = new BufferedReader(new FileReader(new File(fileLocation)));
				String pepsumline = "";
				while((pepsumline=pepsumreader.readLine())!=null) {
					if(pepsumline.startsWith("Peptide Sequence")) {
							String[] pepsumlinesplit = pepsumline.split("\t");
							for(int i = 0; i < pepsumlinesplit.length; ++i) {
								switch(pepsumlinesplit[i]) {
								case "Protein ID":
									prot_acc_idx = i;
									break;
								case "Peptide Sequence":
									pep_seq_idx = i;
									break;
								case "Modified Sequence":
									pep_mod_seq_idx = i;
									break;
								case "Charges":
									pep_exp_z_idx = i;
									break;
								default:
									break;
								}
								if(pepsumlinesplit[i].contains("Intensity")) {
									if(!pepsumlinesplit[i].contains("MaxLFQ Intensity")) {
										exp_idxs.add(i);
										exp_names.add(pepsumlinesplit[i].replaceAll(" Intensity",""));
									}
								}
							}
						} else {
							String[] pepsumlinesplit = pepsumline.split("\t");
							String nterm_cleav = "";
							String cterm_cleav = "";
							
							if(pepsumlinesplit.length>0) {
								if(pepsumlinesplit[prot_acc_idx].contains(id) || (id.equals("P10636") && (pepsumlinesplit[prot_acc_idx].contains("0N3R") || pepsumlinesplit[prot_acc_idx].contains("1N3R") || pepsumlinesplit[prot_acc_idx].contains("2N3R") ||
										pepsumlinesplit[prot_acc_idx].contains("0N4R") || pepsumlinesplit[prot_acc_idx].contains("1N4R") || pepsumlinesplit[prot_acc_idx].contains("2N4R")))) {
									String accession = "";
									int index = -1;
									if(headers_ordered.size()==4) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "2N4R";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
										} else if(taumap.get(headers_ordered.get(1)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "0N";
											index = taumap.get(headers_ordered.get(1)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(1)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(1)),pepsumlinesplit[pep_seq_idx],accession);
										} else if(taumap.get(headers_ordered.get(2)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "1N";
											index = taumap.get(headers_ordered.get(2)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(2)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(2)),pepsumlinesplit[pep_seq_idx],accession);
										} else if(taumap.get(headers_ordered.get(3)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "3R";
											index = taumap.get(headers_ordered.get(3)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(3)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(3)),pepsumlinesplit[pep_seq_idx],accession);
										}
									} else if(headers_ordered.size()==2) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "alpha";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
										} else if(taumap.get(headers_ordered.get(1)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "epsilon";
											index = taumap.get(headers_ordered.get(1)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(1)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(1)),pepsumlinesplit[pep_seq_idx],accession);
										}
									} else if(headers_ordered.size()==3) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "Can";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
										} else if(taumap.get(headers_ordered.get(1)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "Ex5";
											index = taumap.get(headers_ordered.get(1)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(1)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(1)),pepsumlinesplit[pep_seq_idx],accession);
										} else if(taumap.get(headers_ordered.get(2)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "Ex3";
											index = taumap.get(headers_ordered.get(2)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(2)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(2)),pepsumlinesplit[pep_seq_idx],accession);
										}
									} else if(headers_ordered.size()==1) {
										if(taumap.get(headers_ordered.get(0)).contains(pepsumlinesplit[pep_seq_idx])) {
											accession = "PADI2";
											index = taumap.get(headers_ordered.get(0)).indexOf(pepsumlinesplit[pep_seq_idx]);
											nterm_cleav = checkNtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
											cterm_cleav = checkCtermCleavage(taumap.get(headers_ordered.get(0)),pepsumlinesplit[pep_seq_idx],accession);
										}
									}
									String cleavage = "";
									if(!nterm_cleav.equals("") && !cterm_cleav.equals("")) {
										cleavage = nterm_cleav + "; " + cterm_cleav;
									} else if(!nterm_cleav.equals("")) {
										cleavage = nterm_cleav;
									} else if(!cterm_cleav.equals("")) {
										cleavage = cterm_cleav;
									}

									String reformmodsstring = reformmods(pepsumlinesplit[pep_mod_seq_idx],modlookup);
									for(int k = 0; k < exp_idxs.size(); ++k) {
									if(pepsumlinesplit[pep_seq_idx].equals("HVPGGGSVQIVYKPVDLSKVTSK")) {
										boolean test1 = true;
									}
										if(((checkmods(reformmodsstring)==true) || reformmodsstring.equals("")) && !accession.equals("")) {
											writer.write(exp_names.get(k) + sep + accession + sep + "" + sep + pepsumlinesplit[pep_seq_idx] + sep + reformmodsstring + sep + correctingMods(reformmodsstring, index, pepsumlinesplit[pep_seq_idx], accession) + sep);
											writer.write("" + sep + addPTMsequences(pepsumlinesplit[pep_seq_idx], reformmodsstring, reformmodsstring) + sep);
											writer.write("" + sep + cleavage + sep + "" + sep + "" + sep);
											writer.write("" + sep + "" + sep + "" + sep + pepsumlinesplit[pep_exp_z_idx].replaceAll(",",";") + sep + pepsumlinesplit[exp_idxs.get(k)] + "\n");
										}
									}
								}
								
							}
						}
				}
				pepsumreader.close();
//				writer.write(newfile + sep + FDR + sep + threshold + sep + conf + "\n");
	    	}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
    public static void WriteToFile(String str, String filePath) throws IOException {
		BufferedWriter bw = null;
		try {
			 FileOutputStream out = new FileOutputStream(filePath, true); // true, means: file appended content, not regenerated, default is false
			bw = new BufferedWriter(new OutputStreamWriter(out));
			 bw.write(str += "\r\n");// Line break
			bw.flush();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			bw.close();
		}
    }

}
