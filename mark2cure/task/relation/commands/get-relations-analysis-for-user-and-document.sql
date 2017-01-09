SELECT  `relation_relation`.`id` as `relation_id`,
        `relation_relation`.`document_id`,
        `relation_relation`.`relation_type`,

        `relation_relation`.`concept_1_id`,
        (
            SELECT `relation_concepttext`.`text`
            FROM `relation_concepttext`

            INNER JOIN `relation_conceptdocumentrelationship`
              ON `relation_conceptdocumentrelationship`.`document_id` = {document_id}
                AND `relation_conceptdocumentrelationship`.`concept_text_id` = `relation_concepttext`.`id`

            WHERE `relation_concepttext`.`concept_id` = `relation_relation`.`concept_1_id`

            ORDER BY LENGTH(`relation_concepttext`.`text`) DESC
            LIMIT 1
        ) as `concept_1_text`,

        `relation_relation`.`concept_2_id`,
        (
            SELECT `relation_concepttext`.`text`
            FROM `relation_concepttext`

            INNER JOIN `relation_conceptdocumentrelationship`
              ON `relation_conceptdocumentrelationship`.`document_id` = {document_id}
                AND `relation_conceptdocumentrelationship`.`concept_text_id` = `relation_concepttext`.`id`

            WHERE `relation_concepttext`.`concept_id` = `relation_relation`.`concept_2_id`

            ORDER BY LENGTH(`relation_concepttext`.`text`) DESC
            LIMIT 1
        ) as `concept_2_text`,

        `relation_relationannotation`.`answer`,
        `document_view`.`user_id`,
        IF(
            `document_view`.`user_id` = {user_id}
        , TRUE, FALSE) as `self`

FROM `relation_relationannotation`

INNER JOIN `relation_relation`
    ON `relation_relation`.`id` = `relation_relationannotation`.`relation_id`

INNER JOIN `document_annotation`
  ON (`document_annotation`.`object_id` = `relation_relationannotation`.`id`
    AND `document_annotation`.`content_type_id` = 56)

INNER JOIN `document_view`
    ON `document_view`.`id` = `document_annotation`.`view_id`

WHERE `relation_relation`.`document_id` = {document_id}
  AND `relation_id` IN (

      /*  Subquery to select all the relationship identifiers
          for a document that the user has submitted an
          annotation for */
      SELECT  `relation_relation`.`id` as `relation_id`

      FROM `relation_relation`

      INNER JOIN `relation_relationannotation`
          ON `relation_relationannotation`.`relation_id` = `relation_relation`.`id`

      INNER JOIN `document_annotation`
          ON (`document_annotation`.`object_id` = `relation_relationannotation`.`id`
              AND `document_annotation`.`content_type_id` = 56)

      INNER JOIN `document_view`
          ON `document_view`.`id` = `document_annotation`.`view_id`

      WHERE `relation_relation`.`document_id` = {document_id}
          AND `document_view`.`user_id` = {user_id}
    )
    {relation_logic}

ORDER BY  `relation_id` ASC,
          `answer` ASC,
          `user_id` ASC
