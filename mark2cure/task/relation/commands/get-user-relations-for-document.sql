/* Breaks on 5.7.12 */

SELECT  `relationship`.`relation_id`,
        `relationship`.`document_id`,
        `relationship`.`relation_type`,

        `relationship`.`progress`,
        `relationship`.`community_completed`,
        `relationship`.`user_completed`,

        `relationship`.`concept_1_id`,
        (
          SELECT `relation_conceptdocumentrelationship`.`stype`
          FROM `relation_concepttext`

          INNER JOIN `relation_conceptdocumentrelationship`
            ON `relation_conceptdocumentrelationship`.`document_id` = {document_id}
              AND `relation_conceptdocumentrelationship`.`concept_text_id` = `relation_concepttext`.`id`

          WHERE `relation_concepttext`.`concept_id` = `relationship`.`concept_1_id`

          LIMIT 1
        ) as `concept_1_type`,
        (
          SELECT `relation_concepttext`.`text`
          FROM `relation_concepttext`

          INNER JOIN `relation_conceptdocumentrelationship`
            ON `relation_conceptdocumentrelationship`.`document_id` = {document_id}
              AND `relation_conceptdocumentrelationship`.`concept_text_id` = `relation_concepttext`.`id`

          WHERE `relation_concepttext`.`concept_id` = `relationship`.`concept_1_id`

          ORDER BY LENGTH(`relation_concepttext`.`text`) DESC
          LIMIT 1
        ) as `concept_1_text`,

        `relationship`.`concept_2_id`,
        (
          SELECT `relation_conceptdocumentrelationship`.`stype`
          FROM `relation_concepttext`

          INNER JOIN `relation_conceptdocumentrelationship`
            ON `relation_conceptdocumentrelationship`.`document_id` = {document_id}
              AND `relation_conceptdocumentrelationship`.`concept_text_id` = `relation_concepttext`.`id`

          WHERE `relation_concepttext`.`concept_id` = `relationship`.`concept_2_id`

          LIMIT 1
        ) as `concept_2_type`,
        (
          SELECT `relation_concepttext`.`text`
          FROM `relation_concepttext`

          INNER JOIN `relation_conceptdocumentrelationship`
            ON `relation_conceptdocumentrelationship`.`document_id` = {document_id}
              AND `relation_conceptdocumentrelationship`.`concept_text_id` = `relation_concepttext`.`id`

          WHERE `relation_concepttext`.`concept_id` = `relationship`.`concept_2_id`

          ORDER BY LENGTH(`relation_concepttext`.`text`) DESC
          LIMIT 1
        ) as `concept_2_text`

FROM (
  SELECT  `relation_relation`.`id` as `relation_id`,
          `relation_relation`.`document_id`,
          `relation_relation`.`relation_type`,
          `relation_relation`.`concept_1_id`,
          `relation_relation`.`concept_2_id`,

          # (TODO) What happends when there are no annotations?
          # MIN(`document_annotation`.`created`) as `first_annotation`,
          # MAX(`document_annotation`.`created`) as `last_annotation`,

          /*  If the relationship has been answered enough
              times by a minimum number of unique users */
          IF(
            /*  Ensure we count the unique number of users who
                have completed the relationship */
            COUNT(DISTINCT `document_view`.`user_id`) >= {completions}
          , TRUE, FALSE) as `community_completed`,

          /*  If the provided users of interest has provided
              an answer for the specific relationship */
          IF(
              SUM(`document_view`.`user_id` = {user_id})
          , TRUE, FALSE) as `user_completed`,

          ROUND(
            CASE
              WHEN COUNT(DISTINCT `document_view`.`user_id`)/{completions} > 1 THEN 1
              ELSE COUNT(DISTINCT `document_view`.`user_id`)/{completions}
            END
          , 2) as `progress`

  FROM `relation_relation`

  LEFT JOIN `relation_relationannotation`
      ON `relation_relationannotation`.`relation_id` = `relation_relation`.`id`

  LEFT JOIN `document_annotation`
      ON (`document_annotation`.`object_id` = `relation_relationannotation`.`id`
        AND `document_annotation`.`content_type_id` = 56)

  /*  We're only joining Views to get back the User ID,
      We don't care about the "completed" boolean b/c
      the Relation tasks are not all or none */
  LEFT JOIN `document_view`
      ON `document_view`.`id` = `document_annotation`.`view_id`

  WHERE `relation_relation`.`document_id` = {document_id}

  GROUP BY `relation_id`

  ORDER BY `progress` DESC
) as `relationship`

HAVING `relationship`.`community_completed` = FALSE AND NOT `relationship`.`user_completed` = TRUE

LIMIT 20
