-- phpMyAdmin SQL Dump
-- version 4.0.10deb1
-- http://www.phpmyadmin.net
--
-- Host: localhost:3306
-- Generation Time: Jan 19, 2015 at 04:59 PM
-- Server version: 5.5.40-0ubuntu0.14.04.1
-- PHP Version: 5.5.9-1ubuntu4.5

SET FOREIGN_KEY_CHECKS=0;
SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `aspminedb`
--
CREATE DATABASE IF NOT EXISTS `aspminedb` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;
USE `aspminedb`;

-- --------------------------------------------------------

--
-- Table structure for table `antismash`
--
-- Creation: Nov 27, 2014 at 02:22 PM
--

DROP TABLE IF EXISTS `antismash`;
CREATE TABLE IF NOT EXISTS `antismash` (
  `org_id` int(11) NOT NULL,
  `filename` text NOT NULL,
  `clust_size` int(11) NOT NULL,
  `clust_backbone` varchar(50) NOT NULL,
  `clust_origin` varchar(50) NOT NULL,
  `clust_start` int(11) NOT NULL,
  `clust_end` int(11) NOT NULL,
  `clust_type` varchar(50) NOT NULL,
  `sm_id` int(11) NOT NULL,
  `sm_protein_id` int(11) NOT NULL,
  `sm_short` varchar(100) NOT NULL,
  `sm_desc` text NOT NULL,
  `cluster` text NOT NULL,
  PRIMARY KEY (`clust_backbone`,`sm_protein_id`),
  KEY `sm_protein_id` (`sm_protein_id`),
  KEY `clust_backbone` (`clust_backbone`),
  KEY `clust_backbone_sm_protein_id` (`clust_backbone`,`sm_protein_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `assembly`
--
-- Creation: Nov 27, 2014 at 02:22 PM
--

DROP TABLE IF EXISTS `assembly`;
CREATE TABLE IF NOT EXISTS `assembly` (
  `assembly_orgkey` varchar(100) NOT NULL,
  `assembly_seqkey` varchar(100) NOT NULL,
  `assembly_tailkey` varchar(100) NOT NULL,
  `org_id` int(11) NOT NULL,
  `assembly_seqname` varchar(100) NOT NULL,
  `assembly_orgseq` varchar(100) DEFAULT NULL,
  `assembly_seq` longblob,
  PRIMARY KEY (`assembly_seqname`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `blast`
--
-- Creation: Jan 15, 2015 at 09:33 AM
--

DROP TABLE IF EXISTS `blast`;
CREATE TABLE IF NOT EXISTS `blast` (
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_qseq_jg2` varchar(100) NOT NULL,
  `blast_qseq_jg3` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg2` varchar(100) NOT NULL,
  `blast_sseq_jg3` varchar(100) NOT NULL,
  `blast_pident` decimal(10,2) DEFAULT NULL,
  `blast_qlen` int(11) NOT NULL,
  `blast_qstart` int(11) NOT NULL,
  `blast_qend` int(11) NOT NULL,
  `blast_slen` int(11) NOT NULL,
  `blast_sstart` int(11) NOT NULL,
  `blast_send` int(11) NOT NULL,
  `blast_bitscore` int(11) NOT NULL,
  `blast_evalue` varchar(10) NOT NULL,
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_filename` varchar(100) NOT NULL,
  PRIMARY KEY (`blast_bitscore`,`blast_qseq_id`,`blast_sseq_id`,`blast_qstart`),
  KEY `blast_qseq_jg1` (`blast_qseq_jg1`),
  KEY `blast_sseq_jg1` (`blast_sseq_jg1`),
  KEY `blast_qseq_jg3` (`blast_qseq_jg3`),
  KEY `blast_sseq_jg3` (`blast_sseq_jg3`),
  KEY `blast_qseq_jg2` (`blast_qseq_jg2`),
  KEY `blast_sseq_jg2` (`blast_sseq_jg2`),
  KEY `blast_sseq_id` (`blast_sseq_id`),
  KEY `blast_qseq_id` (`blast_qseq_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `blast2org`
--
-- Creation: Nov 27, 2014 at 09:08 PM
--

DROP TABLE IF EXISTS `blast2org`;
CREATE TABLE IF NOT EXISTS `blast2org` (
  `qseq_id` varchar(50) NOT NULL DEFAULT '',
  `sseq_id` varchar(50) NOT NULL DEFAULT '',
  `qorg_id` int(11) DEFAULT NULL,
  `sorg_id` int(11) DEFAULT NULL,
  `qprot_seq_name` varchar(100) DEFAULT NULL,
  `sprot_seq_name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`qseq_id`,`sseq_id`),
  KEY `qorg_id` (`qorg_id`),
  KEY `sorg_id` (`sorg_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `cds`
--
-- Creation: Dec 18, 2014 at 09:16 AM
--

DROP TABLE IF EXISTS `cds`;
CREATE TABLE IF NOT EXISTS `cds` (
  `org_id` int(11) NOT NULL,
  `cds_orgseq` varchar(100) DEFAULT NULL,
  `cds_orgkey` varchar(100) NOT NULL,
  `cds_seqkey` varchar(100) NOT NULL,
  `cds_tailkey` varchar(100) NOT NULL,
  `cds_seqname` varchar(100) NOT NULL,
  `cds_seq` longblob,
  PRIMARY KEY (`cds_seqname`,`org_id`),
  KEY `org_id` (`org_id`),
  KEY `icds_tailkey` (`cds_tailkey`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `environment`
--
-- Creation: Nov 27, 2014 at 09:08 PM
--

DROP TABLE IF EXISTS `environment`;
CREATE TABLE IF NOT EXISTS `environment` (
  `env_id` int(11) NOT NULL AUTO_INCREMENT,
  `env_desc` varchar(100) NOT NULL,
  `env_class` varchar(100) NOT NULL,
  PRIMARY KEY (`env_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `gff`
--
-- Creation: Jan 06, 2015 at 01:55 PM
--

DROP TABLE IF EXISTS `gff`;
CREATE TABLE IF NOT EXISTS `gff` (
  `gff_seqorigin` varchar(100) NOT NULL,
  `org_id` int(11) NOT NULL DEFAULT '0',
  `gff_name` varchar(100) NOT NULL,
  `gff_type` varchar(100) NOT NULL,
  `gff_protein_id` varchar(100) DEFAULT NULL,
  `gff_trans_id` varchar(100) DEFAULT NULL,
  `gff_attributes` varchar(500) NOT NULL DEFAULT '',
  `gff_start` int(11) NOT NULL DEFAULT '0',
  `gff_end` int(11) DEFAULT NULL,
  `gff_length` int(11) DEFAULT NULL,
  `gff_score` varchar(10) DEFAULT NULL,
  `gff_strand` varchar(10) DEFAULT NULL,
  `gff_phase` varchar(10) DEFAULT NULL,
  `gff_exonnr` varchar(500) NOT NULL,
  PRIMARY KEY (`org_id`,`gff_attributes`,`gff_start`),
  KEY `igff_seq_id` (`gff_protein_id`),
  KEY `igff_name` (`gff_name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `go`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `go`;
CREATE TABLE IF NOT EXISTS `go` (
  `go_term_id` varchar(50) NOT NULL DEFAULT '',
  `go_name` varchar(500) DEFAULT NULL,
  `go_termtype` varchar(50) DEFAULT NULL,
  `go_acc` varchar(50) NOT NULL DEFAULT '',
  PRIMARY KEY (`go_term_id`,`go_acc`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `ipr`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `ipr`;
CREATE TABLE IF NOT EXISTS `ipr` (
  `ipr_id` varchar(50) NOT NULL DEFAULT '',
  `ipr_desc` text,
  `ipr_domaindb` varchar(50) DEFAULT NULL,
  `ipr_domain_id` varchar(100) NOT NULL DEFAULT '',
  `ipr_domaindesc` text,
  PRIMARY KEY (`ipr_id`,`ipr_domain_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `kegg`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `kegg`;
CREATE TABLE IF NOT EXISTS `kegg` (
  `kegg_id` int(11) NOT NULL AUTO_INCREMENT,
  `kegg_ecNum` varchar(50) NOT NULL,
  `kegg_definition` varchar(500) DEFAULT NULL,
  `kegg_catalyticActivity` varchar(500) DEFAULT NULL,
  `kegg_cofactors` varchar(500) DEFAULT NULL,
  `kegg_associatedDiseases` varchar(500) DEFAULT NULL,
  `kegg_pathway` varchar(500) NOT NULL,
  `kegg_pathway_class` varchar(500) DEFAULT NULL,
  `kegg_pathway_type` varchar(500) DEFAULT NULL,
  PRIMARY KEY (`kegg_id`),
  UNIQUE KEY `kegg_ecNum_kegg_pathway` (`kegg_ecNum`,`kegg_pathway`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=306963 ;

-- --------------------------------------------------------

--
-- Table structure for table `kog`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `kog`;
CREATE TABLE IF NOT EXISTS `kog` (
  `kog_id` varchar(50) NOT NULL DEFAULT '',
  `kog_defline` varchar(500) DEFAULT NULL,
  `kog_Class` varchar(500) DEFAULT NULL,
  `kog_Group` varchar(500) DEFAULT NULL,
  PRIMARY KEY (`kog_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lengths`
--
-- Creation: Jan 15, 2015 at 03:57 PM
--

DROP TABLE IF EXISTS `lengths`;
CREATE TABLE IF NOT EXISTS `lengths` (
  `org_id` int(11) NOT NULL,
  `prot_seqname` varchar(100) NOT NULL,
  `len` int(10) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `org_has_environment`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `org_has_environment`;
CREATE TABLE IF NOT EXISTS `org_has_environment` (
  `env_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  PRIMARY KEY (`env_id`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `org_is_species`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `org_is_species`;
CREATE TABLE IF NOT EXISTS `org_is_species` (
  `species_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  PRIMARY KEY (`species_id`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `organism`
--
-- Creation: Jan 13, 2015 at 03:25 PM
--

DROP TABLE IF EXISTS `organism`;
CREATE TABLE IF NOT EXISTS `organism` (
  `org_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(100) NOT NULL,
  `real_name` text,
  `section` varchar(100) NOT NULL,
  `status` varchar(50) NOT NULL,
  PRIMARY KEY (`org_id`),
  KEY `iname` (`name`),
  KEY `lable` (`section`)
) ENGINE=InnoDB  DEFAULT CHARSET=latin1 AUTO_INCREMENT=51 ;

-- --------------------------------------------------------

--
-- Table structure for table `primary_metabolism`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `primary_metabolism`;
CREATE TABLE IF NOT EXISTS `primary_metabolism` (
  `react_name` varchar(10) NOT NULL,
  `react_id` varchar(10) NOT NULL,
  `react_path` varchar(100) NOT NULL,
  `reacti_ec` varchar(20) NOT NULL,
  `react_ecname` varchar(100) NOT NULL,
  `react_gene` varchar(20) NOT NULL,
  `react_ref` text NOT NULL,
  `react_cbs51388` text NOT NULL,
  `react_atcc1015` text NOT NULL,
  `react_isoenzymes` varchar(20) NOT NULL,
  `react_modelkey` text NOT NULL,
  PRIMARY KEY (`react_name`,`react_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `primary_metabolism_models`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `primary_metabolism_models`;
CREATE TABLE IF NOT EXISTS `primary_metabolism_models` (
  `model_abb` varchar(100) NOT NULL,
  `model_name` text NOT NULL,
  `model_keggsyn` text NOT NULL,
  `model_keggnr` varchar(10) NOT NULL,
  `model_origin` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_go`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `protein_has_go`;
CREATE TABLE IF NOT EXISTS `protein_has_go` (
  `protein_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  `go_term_id` int(11) NOT NULL,
  PRIMARY KEY (`org_id`,`protein_id`,`go_term_id`),
  KEY `goterm_id` (`go_term_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_ipr`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `protein_has_ipr`;
CREATE TABLE IF NOT EXISTS `protein_has_ipr` (
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  `ipr_id` varchar(50) NOT NULL,
  `ipr_domain_start` int(11) NOT NULL,
  `ipr_domain_end` int(11) NOT NULL,
  `ipr_score` varchar(20) NOT NULL,
  UNIQUE KEY `protein_id_ipr_id` (`protein_id`,`ipr_id`),
  KEY `org_id` (`org_id`),
  KEY `ipr_id` (`ipr_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_kegg`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `protein_has_kegg`;
CREATE TABLE IF NOT EXISTS `protein_has_kegg` (
  `kegg_id` int(11) NOT NULL,
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  `kegg_ecNum` varchar(50) NOT NULL,
  `kegg_pathway` varchar(500) NOT NULL,
  PRIMARY KEY (`kegg_id`,`protein_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `protein_has_kog`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `protein_has_kog`;
CREATE TABLE IF NOT EXISTS `protein_has_kog` (
  `kog_id` varchar(50) NOT NULL,
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  PRIMARY KEY (`kog_id`,`protein_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `proteins`
--
-- Creation: Dec 18, 2014 at 09:13 AM
--

DROP TABLE IF EXISTS `proteins`;
CREATE TABLE IF NOT EXISTS `proteins` (
  `org_id` int(11) NOT NULL,
  `prot_orgkey` varchar(100) NOT NULL,
  `prot_seqkey` varchar(100) NOT NULL,
  `prot_tailkey` varchar(100) NOT NULL,
  `prot_seqname` varchar(100) NOT NULL,
  `prot_orgseq` varchar(100) DEFAULT NULL,
  `prot_seq` longblob,
  PRIMARY KEY (`prot_seqname`,`org_id`),
  KEY `org_id` (`org_id`),
  KEY `iprot_seq_name` (`prot_seqname`),
  KEY `iprot_seqkey` (`prot_seqkey`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `pyrG`
--
-- Creation: Jan 09, 2015 at 10:42 AM
--

DROP TABLE IF EXISTS `pyrG`;
CREATE TABLE IF NOT EXISTS `pyrG` (
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_qseq_jg2` varchar(100) NOT NULL,
  `blast_qseq_jg3` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg2` varchar(100) NOT NULL,
  `blast_sseq_jg3` varchar(100) NOT NULL,
  `blast_pident` decimal(10,2) DEFAULT NULL,
  `blast_qlen` int(11) NOT NULL,
  `blast_qstart` int(11) NOT NULL,
  `blast_qend` int(11) NOT NULL,
  `blast_slen` int(11) NOT NULL,
  `blast_sstart` int(11) NOT NULL,
  `blast_send` int(11) NOT NULL,
  `blast_bitscore` int(11) NOT NULL,
  `blast_evalue` varchar(10) NOT NULL,
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_filename` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `pyrG2`
--
-- Creation: Jan 09, 2015 at 11:30 AM
--

DROP TABLE IF EXISTS `pyrG2`;
CREATE TABLE IF NOT EXISTS `pyrG2` (
  `blast_sseq_id` varchar(100) NOT NULL,
  `idents` text CHARACTER SET utf8,
  `nrTarget` bigint(21) NOT NULL DEFAULT '0',
  `Targets` text,
  `org_id` int(11) NOT NULL,
  `prot_orgkey` varchar(100) NOT NULL,
  `prot_seqkey` varchar(100) NOT NULL,
  `prot_seqname` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `pyrG3`
--
-- Creation: Jan 09, 2015 at 12:53 PM
--

DROP TABLE IF EXISTS `pyrG3`;
CREATE TABLE IF NOT EXISTS `pyrG3` (
  `prot_seqname` varchar(100) NOT NULL,
  `trans_seqkey` varchar(100),
  `trans_seq` longblob
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `rep`
--
-- Creation: Jan 08, 2015 at 10:18 AM
--

DROP TABLE IF EXISTS `rep`;
CREATE TABLE IF NOT EXISTS `rep` (
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_qseq_jg2` varchar(100) NOT NULL,
  `blast_qseq_jg3` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg2` varchar(100) NOT NULL,
  `blast_sseq_jg3` varchar(100) NOT NULL,
  `blast_pident` varchar(50) NOT NULL,
  `blast_qlen` int(11) NOT NULL,
  `blast_qstart` int(11) NOT NULL,
  `blast_qend` int(11) NOT NULL,
  `blast_slen` int(11) NOT NULL,
  `blast_sstart` int(11) NOT NULL,
  `blast_send` int(11) NOT NULL,
  `blast_bitscore` int(11) NOT NULL,
  `blast_evalue` varchar(10) NOT NULL,
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_filename` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `rep2`
--
-- Creation: Jan 08, 2015 at 12:17 PM
--

DROP TABLE IF EXISTS `rep2`;
CREATE TABLE IF NOT EXISTS `rep2` (
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_qseq_jg2` varchar(100) NOT NULL,
  `blast_qseq_jg3` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg2` varchar(100) NOT NULL,
  `blast_sseq_jg3` varchar(100) NOT NULL,
  `blast_pident` varchar(50) NOT NULL,
  `blast_qlen` int(11) NOT NULL,
  `blast_qstart` int(11) NOT NULL,
  `blast_qend` int(11) NOT NULL,
  `blast_slen` int(11) NOT NULL,
  `blast_sstart` int(11) NOT NULL,
  `blast_send` int(11) NOT NULL,
  `blast_bitscore` int(11) NOT NULL,
  `blast_evalue` varchar(10) NOT NULL,
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_filename` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `sigp`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `sigp`;
CREATE TABLE IF NOT EXISTS `sigp` (
  `org_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  `sigp_nn_cutpos` int(11) NOT NULL DEFAULT '0',
  `sigp_neuro_net_vote` int(11) DEFAULT NULL,
  `sigp_hmm_cutpos` int(11) DEFAULT NULL,
  `sigp_hmm_signalpep_probability` float DEFAULT NULL,
  PRIMARY KEY (`org_id`,`protein_id`,`sigp_nn_cutpos`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `species`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `species`;
CREATE TABLE IF NOT EXISTS `species` (
  `species_id` int(11) NOT NULL AUTO_INCREMENT,
  `species_section` int(11) NOT NULL,
  `species_clade` int(11) NOT NULL,
  `species_desc` varchar(100) NOT NULL,
  PRIMARY KEY (`species_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `tmp`
--
-- Creation: Jan 09, 2015 at 02:22 PM
--

DROP TABLE IF EXISTS `tmp`;
CREATE TABLE IF NOT EXISTS `tmp` (
  `blast_qseq_jg1` varchar(100) NOT NULL,
  `blast_qseq_jg2` varchar(100) NOT NULL,
  `blast_qseq_jg3` varchar(100) NOT NULL,
  `blast_sseq_jg1` varchar(100) NOT NULL,
  `blast_sseq_jg2` varchar(100) NOT NULL,
  `blast_sseq_jg3` varchar(100) NOT NULL,
  `blast_pident` decimal(10,2) DEFAULT NULL,
  `blast_qlen` int(11) NOT NULL,
  `blast_qstart` int(11) NOT NULL,
  `blast_qend` int(11) NOT NULL,
  `blast_slen` int(11) NOT NULL,
  `blast_sstart` int(11) NOT NULL,
  `blast_send` int(11) NOT NULL,
  `blast_bitscore` int(11) NOT NULL,
  `blast_evalue` varchar(10) NOT NULL,
  `blast_sseq_id` varchar(100) NOT NULL,
  `blast_qseq_id` varchar(100) NOT NULL,
  `blast_filename` varchar(100) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `tmp2`
--
-- Creation: Jan 09, 2015 at 02:28 PM
--

DROP TABLE IF EXISTS `tmp2`;
CREATE TABLE IF NOT EXISTS `tmp2` (
  `blast_sseq_id` varchar(100) NOT NULL,
  `idents` text CHARACTER SET utf8,
  `nrTarget` bigint(21) NOT NULL DEFAULT '0',
  `Targets` text,
  `org_id` int(11),
  `prot_orgkey` varchar(100),
  `prot_seqkey` varchar(100),
  `prot_seqname` varchar(100)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `tmp3`
--
-- Creation: Jan 09, 2015 at 02:30 PM
--

DROP TABLE IF EXISTS `tmp3`;
CREATE TABLE IF NOT EXISTS `tmp3` (
  `prot_seqname` varchar(100),
  `trans_seqkey` varchar(100),
  `trans_seq` longblob
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `transcripts`
--
-- Creation: Nov 27, 2014 at 09:09 PM
--

DROP TABLE IF EXISTS `transcripts`;
CREATE TABLE IF NOT EXISTS `transcripts` (
  `org_id` int(11) NOT NULL DEFAULT '0',
  `trans_orgseq` longblob,
  `trans_orgkey` varchar(100) NOT NULL,
  `trans_seqkey` varchar(100) NOT NULL,
  `trans_tailkey` varchar(100) NOT NULL,
  `trans_seqname` varchar(100) NOT NULL,
  `trans_seq` longblob,
  PRIMARY KEY (`trans_seqname`,`org_id`),
  KEY `org_id` (`org_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Constraints for dumped tables
--

--
-- Constraints for table `assembly`
--
ALTER TABLE `assembly`
  ADD CONSTRAINT `assembly_ibfk_2` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `blast2org`
--
ALTER TABLE `blast2org`
  ADD CONSTRAINT `blast2org_ibfk_1` FOREIGN KEY (`qorg_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE,
  ADD CONSTRAINT `blast2org_ibfk_2` FOREIGN KEY (`sorg_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

--
-- Constraints for table `cds`
--
ALTER TABLE `cds`
  ADD CONSTRAINT `cds_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

--
-- Constraints for table `gff`
--
ALTER TABLE `gff`
  ADD CONSTRAINT `gff_ibfk_2` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `org_has_environment`
--
ALTER TABLE `org_has_environment`
  ADD CONSTRAINT `org_has_environment_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON UPDATE CASCADE,
  ADD CONSTRAINT `org_has_environment_ibfk_2` FOREIGN KEY (`env_id`) REFERENCES `environment` (`env_id`) ON UPDATE CASCADE;

--
-- Constraints for table `org_is_species`
--
ALTER TABLE `org_is_species`
  ADD CONSTRAINT `org_is_species_ibfk_3` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `org_is_species_ibfk_5` FOREIGN KEY (`species_id`) REFERENCES `species` (`species_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `protein_has_go`
--
ALTER TABLE `protein_has_go`
  ADD CONSTRAINT `protein_has_go_ibfk_2` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `protein_has_ipr`
--
ALTER TABLE `protein_has_ipr`
  ADD CONSTRAINT `protein_has_ipr_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE,
  ADD CONSTRAINT `protein_has_ipr_ibfk_2` FOREIGN KEY (`ipr_id`) REFERENCES `ipr` (`ipr_id`) ON DELETE CASCADE;

--
-- Constraints for table `protein_has_kegg`
--
ALTER TABLE `protein_has_kegg`
  ADD CONSTRAINT `protein_has_kegg_ibfk_4` FOREIGN KEY (`kegg_id`) REFERENCES `kegg` (`kegg_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `protein_has_kegg_ibfk_5` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `protein_has_kog`
--
ALTER TABLE `protein_has_kog`
  ADD CONSTRAINT `protein_has_kog_ibfk_5` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `protein_has_kog_ibfk_6` FOREIGN KEY (`kog_id`) REFERENCES `kog` (`kog_id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints for table `proteins`
--
ALTER TABLE `proteins`
  ADD CONSTRAINT `proteins_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

--
-- Constraints for table `sigp`
--
ALTER TABLE `sigp`
  ADD CONSTRAINT `sigp_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;

--
-- Constraints for table `transcripts`
--
ALTER TABLE `transcripts`
  ADD CONSTRAINT `transcripts_ibfk_1` FOREIGN KEY (`org_id`) REFERENCES `organism` (`org_id`) ON DELETE CASCADE;
SET FOREIGN_KEY_CHECKS=1;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
