SELECT  `relation_relation`.`id` as `relationship_id`,
        `relation_relation`.`document_id`,
        `relation_relation`.`relation_type`,
        `relation_relation`.`concept_1_id`,
        `relation_relation`.`concept_2_id`,

        # ANY_VALUE(`concept_relationship_1`.`stype`) as `concept_1_stype`,
        ANY_VALUE(`concept_text_1`.`text`) as `concept_1_text`,

        # ANY_VALUE(`concept_relationship_2`.`stype`) as `concept_2_stype`,
        ANY_VALUE(`concept_text_2`.`text`) as `concept_2_text`

FROM `relation_relation`

INNER JOIN `relation_conceptdocumentrelationship` as `concept_relationship_1`
    ON `concept_relationship_1`.`document_id` = `relation_relation`.`document_id`

INNER JOIN `relation_concepttext` as `concept_text_1`
    ON `concept_text_1`.`concept_id` = `relation_relation`.`concept_1_id`

INNER JOIN `relation_conceptdocumentrelationship` as `concept_relationship_2`
    ON `concept_relationship_2`.`document_id` = `relation_relation`.`document_id`

INNER JOIN `relation_concepttext` as `concept_text_2`
    ON `concept_text_2`.`concept_id` = `relation_relation`.`concept_2_id`

WHERE `relation_relation`.`document_id` = {document_id} {relation_logic}

GROUP BY `relation_relation`.`id`
